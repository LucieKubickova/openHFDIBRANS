/*---------------------------------------------------------------------------*\
                        _   _ ____ ____ _____ _____ _____ _____ _    _  _____
                       | | | |  __|  _ \_   _|  __ \  __ \  _  \ \  | |/  _  \
  ___  _ __   ___ _ __ | |_| | |_ | | | || | | |_/ / |_/ / |_| |  \ | |  |_|_/
 / _ \| '_ \ / _ \ '_ \|  _  |  _|| | | || | |  __ \  _ ||  _  | \ \| |\___  \
| (_) | |_) |  __/ | | | | | | |  | |/ / | |_| |_/ / | \ \ | | | |\ \ |/ |_|  |
 \___/| .__/ \___|_| |_\_| |_\_|  |___/ \___/\____/|_/  \_|| |_|_| \__|\_____/
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                    with R eynolds A veraged N avier S tokes equations          
-------------------------------------------------------------------------------
License
openHFDIBRANS is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIBRANS. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with RAS turbulence modeling approach and
    wall functions (Kubickova and Isoz, 2023)

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#include "ibInterpolation.H"

using namespace Foam;

//---------------------------------------------------------------------------//
ibInterpolation::ibInterpolation
(
    const fvMesh& mesh,
    const volScalarField& body,
    DynamicList<Tuple2<label,label>>& boundaryCells,
    List<Tuple2<scalar,scalar>>& boundaryDists,
    DynamicList<label>& surfaceCells,
    List<scalar>& surfaceDists,
    labelField& isBoundaryCell
)
:
mesh_(mesh),
body_(body),
surfNorm_
(
    IOobject
    (
        "surfNorm",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedVector("zero", dimless/dimLength, vector::zero)
),
yOrthoi_
(
    IOobject
    (
        "yOrthoi",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("zero", dimless, -1.0)
),
boundaryCells_(boundaryCells),
boundaryDists_(boundaryDists),
surfaceCells_(surfaceCells),
surfaceDists_(surfaceDists),
isBoundaryCell_(isBoundaryCell),
HFDIBDEMDict_
(
    IOobject
    (
        "HFDIBDEMDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
fvSchemes_
(
    IOobject
    (
        "fvSchemes",
        "system",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
)
{
	// read HFDIBDEM dictionary
    boundarySearch_ = HFDIBDEMDict_.lookupOrDefault<word>("boundarySearch", "vertex");
    excludeWalls_ = HFDIBDEMDict_.lookupOrDefault<bool>("excludeWalls", false);
    readSurfNorm_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSurfaceNormal", false);
    aveYOrtho_ = HFDIBDEMDict_.lookupOrDefault<bool>("averageYOrtho", false);
    totalYOrthoAve_ = HFDIBDEMDict_.lookupOrDefault<bool>("totalYOrthoAverage", false);
    intSpan_ = readScalar(HFDIBDEMDict_.lookup("interfaceSpan"));
    thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
    aveCoeff_ = HFDIBDEMDict_.lookupOrDefault<scalar>("averagingCoeff", 1.0);
    nAveYOrtho_ = HFDIBDEMDict_.lookupOrDefault<label>("nAveragingYOrtho", 1.0);
    averageV_ = HFDIBDEMDict_.lookupOrDefault<bool>("averageVolume", false);

    // read fvSchemes
    HFDIBInnerSchemes_ = fvSchemes_.subDict("HFDIBSchemes").subDict("innerSchemes");

    // compute average cell volume
    VAve_ = 0.0;
    if (averageV_)
    {
        forAll(mesh_.V(), i)
        {
            VAve_ += mesh_.V()[i];
        }
        VAve_ /= mesh_.V().size();
    }

    // calculate surface normals
    calculateSurfNorm();

    // prepare label field
    isBoundaryCell_.setSize(mesh_.C().size());
}

//---------------------------------------------------------------------------//
ibInterpolation::~ibInterpolation()
{
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateInterpolationPoints
(
)
{
    // prepare boundary interpolation info list
    intInfoListBoundary_.setSize(boundaryCells_.size());

    // loop over boundary cells
    forAll(boundaryCells_, bCell)
    {
        // get origin cell label
        label cellI = boundaryCells_[bCell].first();

        // find surf point
        point surfPoint = mesh_.C()[cellI];
        scalar ds = boundaryDists_[bCell].second();
        vector surfNormToSend = surfNorm_[cellI];
        surfPoint -= surfNormToSend*ds;

        // get interpolation point
        getInterpolationPoint(cellI, surfPoint, surfNormToSend, intInfoListBoundary_[bCell]);
    }

    // set interpolation order
    setInterpolationOrder(intInfoListBoundary_);

    // prepare surface interpolation info list
    intInfoListSurface_.setSize(surfaceCells_.size());

    // loop over surface cells
    forAll(surfaceCells_, sCell)
    {
        // get origin cell label
        label cellI = surfaceCells_[sCell];

        // find surf point
        point surfPoint = mesh_.C()[cellI];
        scalar ds = surfaceDists_[sCell];
        vector surfNormToSend = surfNorm_[cellI];
        surfPoint -= surfNormToSend*ds;

        // get interpolation point
        getInterpolationPoint(cellI, surfPoint, surfNormToSend, intInfoListSurface_[sCell]);
    }

    // set interpolation order
    setInterpolationOrder(intInfoListSurface_);
}

//---------------------------------------------------------------------------//
void ibInterpolation::getInterpolationPoint
(
    label cellI,
    point surfPoint,
    vector surfNormToSend,
    interpolationInfo& intInfo
)
{
   // create vector for points and cells and add to main vectors
   DynamicList<point> intPoints;
   DynamicList<label> intCells;

   // prepare reference distance
   scalar intDist;
   if (averageV_)
   {
       intDist = Foam::pow(VAve_,0.333);
   }
   else
   {
       intDist = Foam::pow(mesh_.V()[cellI],0.333);
   }
   intDist *= 0.5;

   // add to list
   intPoints.append(surfPoint);
   if (mag(surfNormToSend) > SMALL)
   {
       Tuple2<label,label> helpTup(cellI, -1);
       Tuple2<vector,Tuple2<label,label>> startCell(mesh_.C()[cellI],helpTup);

       // add other interpolation points
       for (int order=0;order<ORDER;order++)
       {
           startCell = findCellCustom(startCell.first(),startCell.second().first(),startCell.second().second(),surfNormToSend,intDist);
           surfPoint = startCell.first();
           cellI = startCell.second().first();

           if (startCell.second().second() == -1)
           {
               if (startCell.second().first() != -1)
               {
                   if (body_[cellI] >= 0.5)
                   {
                       order--;
                       continue;
                   }
               }

               intPoints.append(surfPoint);
               intCells.append(cellI);
           }

           else
           {
               intPoints.append(surfPoint);
               intCells.append(-1);
           }
       }
   }

   else
   {
       for (int order=0;order<ORDER;order++)
       {
           intPoints.append(surfPoint);
           intCells.append(-1);
       }
   }
   
   // assign to global variables
   intInfo.intPoints_ = intPoints;
   intInfo.intCells_ = intCells;
}
    
//---------------------------------------------------------------------------//
void ibInterpolation::setInterpolationOrder
(
    List<interpolationInfo>& intInfoList
)
{
    // decide which order should be used
    for (label infoI = 0; infoI < intInfoList.size(); infoI++)
    {
        List<bool> allowedOrder;
        allowedOrder.setSize(ORDER);

        for (int intPoint=0;intPoint<ORDER;intPoint++)
        {
            if (intInfoList[infoI].intCells_[intPoint] == -1)
            {
                allowedOrder[intPoint] = false;
            }

            else
            {
                allowedOrder[intPoint] = true;
            }
        }

        intInfoList[infoI].order_ = 2;
        if ( allowedOrder[1] == false)
        {
            intInfoList[infoI].order_ = 1;
        }

        // check if first order is possible
        if ( allowedOrder[0] == false)
        {
            intInfoList[infoI].order_ = 0;
        }

        if (intInfoList[infoI].order_ == 2)
        {
            if (intInfoList[infoI].intCells_[0] == intInfoList[infoI].intCells_[1])
                intInfoList[infoI].order_ = 1;
        }
    }
}

//---------------------------------------------------------------------------//
// Custom function to find cell containing point
// Note: this is much cheaper than standard OF functions
Tuple2<vector,Tuple2<label,label>> ibInterpolation::findCellCustom
(
    vector& prevPoint,
    label& startCell,
    label& startProc,
    vector& gradToBody,
    scalar& intDist
)
{
    if (startProc != -1)
    {

        Tuple2<label,label> helpTup(startCell,startProc);
        Tuple2<vector,Tuple2<label,label>> tupleToReturn(prevPoint + intDist*gradToBody,helpTup);
        return tupleToReturn;
    }

    else if (startCell == -1)
    {
        Tuple2<label,label> helpTup(-1,-1);
        Tuple2<vector,Tuple2<label,label>> tupleToReturn(prevPoint + intDist*gradToBody,helpTup);
        return tupleToReturn;
    }

    labelList cellFaces(mesh_.cells()[startCell]);
    label bestFace(0);

    scalar dotProd(-GREAT);
    forAll (cellFaces,faceI)
    {
        vector vecI(mesh_.Cf()[cellFaces[faceI]] - mesh_.C()[startCell]);
        vecI /= mag(vecI);
        scalar auxDotProd(vecI & gradToBody);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            bestFace = cellFaces[faceI];
        }
    }

    labelList cellPoints(mesh_.faces()[bestFace]);
    DynamicList<Tuple2<vector,Tuple2<label,label>>> pointsToCheck;

    forAll (cellPoints, pointI)
    {
        labelList pointFaces(mesh_.pointFaces()[cellPoints[pointI]]);

        forAll (pointFaces, faceI)
        {
            if (mesh_.isInternalFace(pointFaces[faceI]))
            {
                label owner(mesh_.owner()[pointFaces[faceI]]);
                label neighbour(mesh_.neighbour()[pointFaces[faceI]]);

                if (owner != startCell)
                {
                    bool add(true);

                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].second().first() == owner)
                        {
                            add = false;
                            break;
                        }
                    }

                    if (add)
                    {
                        Tuple2<label,label> helpTup(owner,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(mesh_.C()[owner],helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }

                if (neighbour != startCell)
                {
                    bool add(true);
                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].second().first() == neighbour)
                        {
                            add = false;
                            break;
                        }
                    }

                    if (add)
                    {
                        Tuple2<label,label> helpTup(neighbour,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(mesh_.C()[neighbour],helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
            }

            else
            {
                label owner(mesh_.faceOwner()[pointFaces[faceI]]);
                vector distToFace(mesh_.Cf()[pointFaces[faceI]] - mesh_.C()[owner]);
                vector sfUnit(mesh_.Sf()[pointFaces[faceI]]/mag(mesh_.Sf()[pointFaces[faceI]]));
                vector pointToAppend(mesh_.C()[owner] + mag(distToFace)*sfUnit);

                bool add(true);
                forAll (pointsToCheck, i)
                {
                    if (pointsToCheck[i].first() == pointToAppend)
                    {
                        add = false;
                        break;
                    }
                }

                if (add)
                {
                    Tuple2<label,label> helpTup(-1,-1);
                    Tuple2<vector,Tuple2<label,label>> tupleToAdd(pointToAppend,helpTup);
                    pointsToCheck.append(tupleToAdd);
                }
            }
        }
    }

    dotProd = -GREAT;
    Tuple2<label,label> helpTup(-1,-1);
    Tuple2<vector,Tuple2<label,label>> tupleToReturn(vector::zero,helpTup);

    forAll (pointsToCheck,pointI)
    {
        vector vecI(pointsToCheck[pointI].first() - mesh_.C()[startCell]);
        vecI /= (mag(vecI)+SMALL);
        scalar auxDotProd(vecI & gradToBody);

        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            tupleToReturn   = pointsToCheck[pointI];
        }
    }

    if (tupleToReturn.second().first() != -1 and tupleToReturn.second().second() == -1)
    {
        // create interpolation point 
        intDist = mag(prevPoint - tupleToReturn.first());
        point iP = prevPoint + intDist*gradToBody;
        tupleToReturn.first() = iP;

        label cellI = tupleToReturn.second().first();

        // create list of all cells to check
        List<label> cellIs(1 + mesh_.cellCells()[cellI].size());

        cellIs[0] = cellI;
        forAll(mesh_.cellCells()[cellI], nI)
        {
            cellIs[nI+1] = mesh_.cellCells()[cellI][nI];
        }

        // check whether the point is in the already found cell or its neighbours
        forAll(cellIs, cI)
        {
            cellI = cellIs[cI];

            bool isInside(true);

            forAll(mesh_.cells()[cellI], fI)
            {
                label faceI = mesh_.cells()[cellI][fI];
                label owner = mesh_.faceOwner()[faceI];

                vector outNorm = mesh_.Sf()[faceI];

                if (owner != cellI)
                {
                    outNorm *= -1;
                }

                scalar auxDot = (iP - mesh_.Cf()[faceI]) & outNorm;

                if (auxDot > 0)
                {
                    isInside = false;
                }
            }

            if (isInside)
            {
                tupleToReturn.second().first() = cellI;
                break;
            }
        }
    }
    return tupleToReturn;
}

//---------------------------------------------------------------------------//
void ibInterpolation::findBoundaryCells
(
)
{
    // get wall patches
    DynamicList<label> wPatchIs;
    forAll(mesh_.boundary(), pI)
    {
        if (mesh_.boundary()[pI].type() == "wall")
        {
            wPatchIs.append(pI);
        }
    }

    // preparation
    Tuple2<label,label> toAppend;
    isBoundaryCell_ = -1;

    // loop over cells
    forAll(mesh_.cellCells(), cellI)
    {
        bool toInclude(false);
        label vertex(-1);

        if (body_[cellI] < 0.5 && body_[cellI] >= thrSurf_)
        {
            toInclude = true;
        }

        else if (body_[cellI] < thrSurf_)
        {
            if (boundarySearch_ == "vertex")
            {
                forAll(mesh_.cellPoints()[cellI], pID) // vertex neighbours
                {
                    label pointI = mesh_.cellPoints()[cellI][pID];

                    forAll(mesh_.pointCells()[pointI], cI)
                    {
                        if (body_[mesh_.pointCells()[pointI][cI]] >= 0.5)
                        {
                            toInclude = true;
                            vertex = mesh_.pointCells()[pointI][cI];
                            break;
                        }
                    }
                }
            }

            else if (boundarySearch_ == "face")
            {
                forAll(mesh_.cellCells()[cellI], nID) // face neighbours
                {
                    if (body_[mesh_.cellCells()[cellI][nID]] >= 0.5)
                    {
                        toInclude = true;
                        vertex = mesh_.cellCells()[cellI][nID];
                        break;
                    }
                }
            }

            else
            {
                FatalError << "Boundary cells search " << boundarySearch_ << " not implemented" << exit(FatalError);
            }
        }

        // check whether the cell is adjecent to a regular wall
        if (toInclude and excludeWalls_)
        {
            forAll(mesh_.cells()[cellI], f)
            {
                // get face label
                label faceI = mesh_.cells()[cellI][f];

                if (faceI >= mesh_.owner().size())
                {
                    bool wallFace(false);

                    // loop over patches of type wall
                    forAll(wPatchIs, pI)
                    {
                        // get patch label
                        label patchI = wPatchIs[pI];

                        // get start and end face index
                        label startI = mesh_.boundary()[patchI].start();
                        label endI = startI + mesh_.boundary()[patchI].Cf().size();

                        if (faceI >= startI and faceI < endI)
                        {
                            wallFace = true;
                        }
                    }

                    // exclude wall faces
                    if (wallFace)
                    {
                        toInclude = false;
                    }
                }
            }
        }

        // add the cell
        if (toInclude)
        {
            if (vertex == -1)
            {
                Tuple2<label,label> helpTup(cellI,-1);
                Tuple2<vector,Tuple2<label,label>> startCell(mesh_.C()[cellI],helpTup);

                scalar intDist;
                if (averageV_)
                {
                    intDist = Foam::pow(VAve_,0.333);
                }
                else
                {
                    intDist = Foam::pow(mesh_.V()[cellI],0.333);
                }
                intDist *= 0.5;

                vector surfNormToSend(-surfNorm_[cellI]);
                startCell = findCellCustom(startCell.first(),startCell.second().first(),startCell.second().second(),surfNormToSend,intDist);

                toAppend.first() = cellI;
                toAppend.second() = startCell.second().first();
            }

            else
            {
                toAppend.first() = cellI;
                toAppend.second() = vertex;
            }

            // fix non-existent inner cells
            if (toAppend.second() == -1)
            {
                toAppend.second() = toAppend.first();
            }

            boundaryCells_.append(toAppend);
        }
    }

    // prepare label field
    forAll(boundaryCells_, bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[bCell].first();

        // assign
        isBoundaryCell_[cellI] = bCell;
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::findSurfaceCells
(
)
{
    // loop over all cells
    forAll(body_, cellI)
    {
        // check lambda
        if (body_[cellI] >= thrSurf_)
        {
            if (body_[cellI] <= (1 - thrSurf_))
            {
                surfaceCells_.append(cellI);
            }
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::setUpSurface
(
    volScalarField& surface,
    scalar boundaryValue
)
{
    // reset field
    surface *= 0.0;

    // find in-solid cells
    forAll(body_, cellI)
    {
        if (body_[cellI] >= 0.5)
        {
            surface[cellI] = 1.0;
        }
    }

    // find boundary cells
    forAll(boundaryCells_, bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[bCell].first();

        // set the surface value
        surface[cellI] = boundaryValue;
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::setLambdaBasedSurface
(
    volScalarField& surface,
    scalar boundaryValue
)
{
    // reset field
    surface *= 0.0;

    // loop over cells
    forAll(body_, cellI)
    {
        if (body_[cellI] >= boundaryValue)
        {
            surface[cellI] = 1.0;
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::updateSwitchSurface
(
    volScalarField& surface,
    volScalarField& yPlusi,
    scalar yPlusLam
)
{
    // update surface based on yPlus value
    forAll(boundaryCells_, bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[bCell].first();

        // set the surface value
        if (yPlusi[cellI] <= yPlusLam)
        {
            surface[cellI] = 1.0;
        }

        else
        {
            surface[cellI] = 0.0;
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateSurfNorm
(
)
{
    if (not readSurfNorm_)
    {
        // stabilisation
        dimensionedScalar deltaN("deltaN", dimless/dimLength, SMALL);

        // calculate the surface normal based on the body gradient
        surfNorm_ = -fvc::grad(body_);
        surfNorm_ /= (mag(surfNorm_) + deltaN);
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateBoundaryDist
(
)
{
    // loop over boundary cells
    forAll(boundaryCells_, bCell)
    {
        // get the outer and innter cell label
        label outCellI = boundaryCells_[bCell].first();
        label inCellI = boundaryCells_[bCell].second();

        // prepare
        Tuple2<scalar,scalar> toSave;
        scalar sigma;
        scalar yOrtho;
        point surfPoint;
        scalar l;

        // if outer cell is intersected
        if (body_[outCellI] >= thrSurf_)
        {
            if (averageV_)
            {
                l = Foam::pow(VAve_, 0.333);
            }
            else
            {
                l = Foam::pow(mesh_.V()[outCellI], 0.333);
            }
            sigma = Foam::atanh(1-2*body_[outCellI])*l/intSpan_; // y > 1 for lambda < 0.5
            //~ yOrtho = sigma; // standard approach
            yOrtho = 0.5*(sigma + l*0.5);
        }

        // if inner cell is intersected
        else if (body_[inCellI] < 1.0)
        {
            if (averageV_)
            {
                l = Foam::pow(VAve_, 0.333);
            }
            else
            {
                l = Foam::pow(mesh_.V()[inCellI], 0.333);
            }
            sigma = -1*Foam::atanh(1-2*body_[inCellI])*l/intSpan_; // y > 1 for lambda < 0.5
            surfPoint = mesh_.C()[inCellI];
            surfPoint += surfNorm_[inCellI]*sigma;
            //~ yOrtho = sigma; // standard approach
            yOrtho = surfNorm_[inCellI] & (mesh_.C()[outCellI] - surfPoint);
            yOrtho = 0.5*(yOrtho + l*0.5);
        }

        // not intersected outer cells
        else
        {
            // find the shared vertices
            DynamicList<label> sharedVers;
            forAll(mesh_.cellPoints()[inCellI], iI)
            {
                forAll(mesh_.cellPoints()[outCellI], oI)
                {
                    if (mesh_.cellPoints()[inCellI][iI] == mesh_.cellPoints()[outCellI][oI])
                    {
                        sharedVers.append(mesh_.cellPoints()[inCellI][iI]);
                    }
                }
            }

            // average vertices
            vector center(vector::zero);
            forAll(sharedVers, vI)
            {
                center += mesh_.points()[sharedVers[vI]];
            }
            center /= sharedVers.size();

            // compute ds as a distance from the cell center to the averaged vertex
            sigma = mag(mesh_.C()[inCellI] - center);
            yOrtho = mag(mesh_.C()[outCellI] - center);
        }

        toSave.first() = yOrtho;
        toSave.second() = sigma;
        boundaryDists_[bCell] = toSave;
    }

    // average orthogonal distance over all boundary cells
    if (totalYOrthoAve_)
    {
        // preparation
        scalar yAve(0.0);

        // loop over boundary cells
        forAll(boundaryCells_, bCell)
        {
            yAve += boundaryDists_[bCell].first();
        }

        // calculate average
        yAve /= boundaryDists_.size();

        // save to all boundary cells
        forAll(boundaryCells_, bCell)
        {
            boundaryDists_[bCell].first() = yAve;
        }
    }

    // average orthogonal distance over neighbors
    else if (aveYOrtho_)
    {
        for (label nCorr = 0; nCorr < nAveYOrtho_; nCorr++)
        {
            // save old values
            List<Tuple2<scalar,scalar>> yOld = boundaryDists_;

            // loop over boundary cells
            forAll(boundaryCells_, bCell)
            {
                // get the cell label
                label cellI = boundaryCells_[bCell].first();

                // prepare
                scalar yAve(0.0);
                scalar scalesSum(0.0);
                label nAdded(0);

                // loop over vertices
                forAll(mesh_.cellPoints()[cellI], pID)
                {
                    // get the point label
                    label pointI = mesh_.cellPoints()[cellI][pID];

                    // loop over vertex neighbors
                    forAll(mesh_.pointCells()[pointI], cI)
                    {
                        // get the neighbor cell label
                        label cellN = mesh_.pointCells()[pointI][cI];
                        
                        // skip current cell
                        if (cellI == cellN)
                        {
                            continue;
                        }

                        // skip non boundary cells
                        if (isBoundaryCell_[cellN] == -1)
                        {
                            continue;
                        }

                        // get local boundary cell label
                        label nCell = isBoundaryCell_[cellN];

                        // calculate dot product of surf normals
                        scalar dotNorm = surfNorm_[cellI] & surfNorm_[cellN];

                        // skip if dotNorm negative
                        if (dotNorm < 0.0)
                        {
                            continue;
                        }

                        // calculate distance of cell centers
                        scalar delta = mag(mesh_.C()[cellI] - mesh_.C()[cellN]);

                        // add to average y and scales sum
                        yAve += dotNorm/delta*boundaryDists_[nCell].first();
                        scalesSum += dotNorm/delta;
                        nAdded += 1;
                    }
                }

                // compute average scale
                scalar aveScale = scalesSum/nAdded;

                // add current cell
                yAve += aveCoeff_*aveScale*boundaryDists_[bCell].first();
                scalesSum += aveCoeff_*aveScale;

                // divide average y by scales sum
                yAve /= scalesSum;

                // assign
                boundaryDists_[bCell].first() = yAve;
            }
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateSurfaceDist
(
)
{
    // prepare
    scalar sigma;

    // loop over surface cells
    forAll(surfaceCells_, sCell)
    {
        // get cell label
        label cellI = surfaceCells_[sCell];

        // get cell dimension
        scalar l;
        if (averageV_)
        {
            l = Foam::pow(VAve_, 0.333);
        }
        else
        {
            l = Foam::pow(mesh_.V()[cellI], 0.333);
        }

        // calculate signed distance
        sigma = Foam::atanh(1-2*body_[cellI])*l/intSpan_;

        // assign
        surfaceDists_[sCell] = sigma;
        yOrthoi_[cellI] = sigma;
    }
}

//---------------------------------------------------------------------------//
autoPtr<ibScheme> ibInterpolation::chosenInterpFunc
(
    word name
)
{
    autoPtr<ibScheme> funcPtr;

    if (name == "constant")
    {
        funcPtr.set(new constantScheme());
    }

    else if (name == "linear")
    {
        funcPtr.set(new linearScheme());
    }

    else if (name == "quadratic")
    {
        funcPtr.set(new quadraticScheme());
    }

    else if (name == "logarithmic")
    {
        funcPtr.set(new logarithmicScheme());
    }

    else if (name == "fixedGradient")
    {
        funcPtr.set(new fixedGradientScheme());
    }

    else if (name == "zeroGradient")
    {
        funcPtr.set(new zeroGradientScheme());
    }

    else
    {
        FatalError << "Interpolation function " << name << " not implemented" << exit(FatalError);
    }

    return funcPtr;
}

//---------------------------------------------------------------------------//
void ibInterpolation::saveInterpolationInfo
(
    word outDir,
    word name
)
{
    // prepare file
    autoPtr<OFstream> outFilePtr;
    word fileName = name + "_boundary.dat";
    outFilePtr.reset(new OFstream(outDir/fileName));
    outFilePtr() << "cellI,inCellI,cellCenter,surfNorm,surfNormIn,bodyOut,bodyIn,yOrtho,sigma,order,intPoints,intCells" << endl;

    // loop over boundary cells
    forAll(boundaryCells_, bCell)
    {
        outFilePtr() << boundaryCells_[bCell].first() << ","
            << boundaryCells_[bCell].second() << ","
            << mesh_.C()[boundaryCells_[bCell].first()] << ","
            << surfNorm_[boundaryCells_[bCell].first()] << ","
            << surfNorm_[boundaryCells_[bCell].second()] << ","
            << body_[boundaryCells_[bCell].first()] << ","
            << body_[boundaryCells_[bCell].second()] << ","
            << boundaryDists_[bCell].first() << ","
            << boundaryDists_[bCell].second() << ","
            << intInfoListBoundary_[bCell].order_ << ","
            << intInfoListBoundary_[bCell].intPoints_ << ","
            << intInfoListBoundary_[bCell].intCells_ << endl;
    }

    // prepare file
    fileName = name + "_surface.dat";
    outFilePtr.reset(new OFstream(outDir/fileName));
    outFilePtr() << "cellI,cellCenter,surfNorm,body,sigma,order,intPoints,intCells" << endl;

    // loop over surface cells
    forAll(surfaceCells_, sCell)
    {
        outFilePtr () << surfaceCells_[sCell] << ","
            << mesh_.C()[surfaceCells_[sCell]] << ","
            << surfNorm_[surfaceCells_[sCell]] << ","
            << body_[surfaceCells_[sCell]] << ","
            << surfaceDists_[sCell] << ","
            << intInfoListSurface_[sCell].order_ << ","
            << intInfoListSurface_[sCell].intPoints_ << ","
            << intInfoListSurface_[sCell].intCells_ << endl;
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::saveBoundaryCells
(
)
{
    List<label> saveOutCells(boundaryCells_.size());
    List<label> saveInCells(boundaryCells_.size());

    forAll(boundaryCells_, bCell)
    {
        saveOutCells[bCell] = boundaryCells_[bCell].first();
        saveInCells[bCell] = boundaryCells_[bCell].second();
    }

    saveCellSet(saveOutCells, "outerBoundaryCells");
    saveCellSet(saveInCells, "innerBoundaryCells");
}

//---------------------------------------------------------------------------//
void ibInterpolation::saveSurfaceCells
(
)
{
    List<label> saveCells(surfaceCells_.size());

    forAll(surfaceCells_, sCell)
    {
        saveCells[sCell] = surfaceCells_[sCell];
    }

    saveCellSet(saveCells, "surfaceCells");
}

//---------------------------------------------------------------------------//
void ibInterpolation::saveCellSet
(
    List<label>& listToSave,
    word fileName
)
{
    cellSet setToSave
    (
        IOobject
        (
            fileName,
            polyMesh::meshSubDir/"sets",
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    word outDir = "./constant/polyMesh/sets";
    if (!isDir(outDir))
    {
        mkDir(outDir);
    }

    OFstream outFile(outDir + "/" + fileName);
    setToSave.writeHeader(outFile);
    outFile << listToSave << endl;
    outFile << "// ************************************************************************* //";
}

//---------------------------------------------------------------------------//
void ibInterpolation::cutFInBoundaryCells
(
    volVectorField& f
)
{
    // loop over boundary cells
    forAll(boundaryCells_, bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[bCell].first();

        // cut f
        f[cellI] = (f[cellI] & surfNorm_[cellI])*surfNorm_[cellI];
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::cutUInBoundaryCells
(
    volVectorField& U
)
{
    // loop over boundary cells
    forAll(boundaryCells_, bCell)
    {
        // get the cell labels
        label cellI = boundaryCells_[bCell].first();
        //~ label outCellI = boundaryCells_[bCell].first();
        //~ label inCellI = boundaryCells_[bCell].second();

        // compute scalar product
        scalar prodUNorm = U[cellI] & surfNorm_[cellI];
        //~ scalar prodUNorm = U[inCellI] & surfNorm_[outCellI];

        // cut U if it aims in opposite to surfNorm
        if (prodUNorm < 0.0)
        {
            U[cellI] -= (U[cellI] & surfNorm_[cellI])*surfNorm_[cellI];
            //~ U[inCellI] -= (U[inCellI] & surfNorm_[outCellI])*surfNorm_[outCellI];
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::cutPhiInBoundaryCells
(
    surfaceScalarField& phi
)
{
    // loop over cells
    forAll(mesh_.C(), cellI)
    {
        // skip outer cells
        if (body_[cellI] < 0.5)
        {
            continue;
        }
        
        // check if it is a boundary cell
        bool isBoundary = false;
        forAll(boundaryCells_, bCell)
        {
            // get the cell label
            label cellB = boundaryCells_[bCell].second();

            if (cellI == cellB)
            {
                isBoundary = true;
            }
        }

        // loop over faces
        forAll(mesh_.cells()[cellI], fI)
        {
            // get the face label
            label faceI = mesh_.cells()[cellI][fI];

            // skip non-internal faces
            if (faceI > mesh_.nInternalFaces())
            {
                continue;
            }

            if (isBoundary)
            {
                // get the face normal
                vector Sf = mesh_.Sf()[faceI];

                // turn off phi for inner faces
                if ((Sf & surfNorm_[cellI]) < 0.0)
                {
                    phi[faceI] = 1e-20;
                }
            }

            else
            {
                phi[faceI] = 1e-20;
            }
        }
    }
}

// ************************************************************************* //
