/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#include "ibInterpolation.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
ibInterpolation::ibInterpolation
(
    const fvMesh& mesh,
    const volScalarField& body,
    DynamicList<Tuple2<label,label>>& boundaryCells,
    List<Tuple2<scalar,scalar>>& boundaryDists,
    DynamicList<label>& boundaryFaces,
    List<bool>& isWallCell
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
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedVector("zero", dimless/dimLength, vector::zero)
),
surfTan_
(
    IOobject
    (
        "surfTan",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedVector("zero", dimless, vector::zero)
),
boundaryCells_(boundaryCells),
boundaryDists_(boundaryDists),
boundaryFaces_(boundaryFaces),
isWallCell_(isWallCell),
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
)
{
	// read HFDIBDEM dictionary
    HFDIBInterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");
    intSpan_ = readScalar(HFDIBDEMDict_.lookup("interfaceSpan"));
    thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));

    // compute average cell volume
    VAve_ = 0.0;
    forAll(mesh_.V(), i)
    {
        VAve_ += mesh_.V()[i];
    }
    VAve_ /= mesh_.V().size();

    // calculate surface normals
    calculateSurfNorm();
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
    // prepare interpolation info list
    intInfoList_.setSize(boundaryCells_.size());

    forAll(boundaryCells_, bCell)
    {
        // get origin cell label
        label cellI = boundaryCells_[bCell].first();

        // find surf point
        point surfPoint = mesh_.C()[cellI];
        scalar ds = boundaryDists_[bCell].first();
        vector surfNormToSend = surfNorm_[cellI];
        surfPoint -= surfNormToSend*ds;

        // create vector for points and cells and add to main vectors
        DynamicList<point> intPoints;
        DynamicList<label> intCells;

        scalar intDist = Foam::pow(VAve_,0.333);
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
        intInfoList_[bCell].intPoints_ = intPoints;
        intInfoList_[bCell].intCells_ = intCells;
    }
    
    // decide which order should be used
    for (label infoI = 0; infoI < intInfoList_.size(); infoI++)
    {
        List<bool> allowedOrder;
        allowedOrder.setSize(ORDER);

        for (int intPoint=0;intPoint<ORDER;intPoint++)
        {
            if (intInfoList_[infoI].intCells_[intPoint] == -1)
            {
                allowedOrder[intPoint] = false;
            }

            else
            {
                allowedOrder[intPoint] = true;
            }
        }

        intInfoList_[infoI].order_ = 2;
        if ( allowedOrder[1] == false)
        {
            intInfoList_[infoI].order_ = 1;
        }

        // check if first order is possible
        if ( allowedOrder[0] == false)
        {
            intInfoList_[infoI].order_ = 0;
        }

        if (intInfoList_[infoI].order_ == 2)
        {
            if (intInfoList_[infoI].intCells_[0] == intInfoList_[infoI].intCells_[1])
                intInfoList_[infoI].order_ = 1;
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
    // get label of wallInsideLambda patch and starting and ending face index
    label patchIL = mesh_.boundaryMesh().findPatchID("wallInsideLambda");
    label startIL = mesh_.boundary()[patchIL].start();
    label endIL = startIL + mesh_.boundary()[patchIL].Cf().size();

    label patchIW = mesh_.boundaryMesh().findPatchID("walls");
    label startIW = mesh_.boundary()[patchIW].start();
    label endIW = startIW + mesh_.boundary()[patchIW].Cf().size();

    // preparation
    Tuple2<label,label> toAppend;

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

            //~ forAll(mesh_.cellCells()[cellI], nID) // face neighbours
            //~ {
                //~ if (body_[mesh_.cellCells()[cellI][nID]] >= 0.5)
                //~ {
                    //~ toInclude = true;
                    //~ break;
                //~ }
            //~ }
        }

        // check whether the cell is adjecent to a regular wall
        if (toInclude)
        {
            forAll(mesh_.cells()[cellI], f)
            {
                // get face label
                label faceI = mesh_.cells()[cellI][f];

                if ((faceI >= startIL and faceI < endIL) or (faceI >= startIW and faceI < endIW))
                {
                    toInclude = false;
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

                scalar intDist = Foam::pow(VAve_,0.333);
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

            boundaryCells_.append(toAppend);
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::areWallCells
(
)
{
    // get label of wallInsideLambda patch and starting and ending face index
    label patchI = mesh_.boundaryMesh().findPatchID("wallInsideLambda");
    label startI = mesh_.boundary()[patchI].start();
    label endI = startI + mesh_.boundary()[patchI].Cf().size();

    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // initialize value
        isWallCell_[bCell] = false;

        // check whether the cell is adjecent to a regular wall
        forAll(mesh_.cells()[cellI], f)
        {
            // get face label
            label faceI = mesh_.cells()[cellI][f];

            if (faceI >= startI and faceI < endI)
            {
                isWallCell_[bCell] = true;
                break;
            }
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::findBoundaryFaces
(
)
{
    forAll(mesh_.cellCells(), cellI)
    {
        label sharedFace = -1;

        if (body_[cellI] >= 0.5)
        {
            forAll(mesh_.cellCells()[cellI], nCell)
            {
                // get the cell label of the neighbor
                label nI = mesh_.cellCells()[cellI][nCell];

                if (body_[nI] < 0.5)
                {
                    // find the shared face
                    forAll(mesh_.cells()[cellI], iFace)
                    {
                        forAll(mesh_.cells()[nI], oFace)
                        {
                            if (iFace == oFace)
                            {
                                sharedFace = mesh_.cells()[cellI][iFace];
                                break;

                            }
                        }

                        if (sharedFace > -1)
                        {
                            break;
                        }
                    }
                }

                if (sharedFace > -1)
                {
                    break;
                }
            }
        }

        if (sharedFace > -1)
        {
            boundaryFaces_.append(sharedFace);
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::createOuterSurface
(
    volScalarField& outSurface
)
{
    // reset field
    outSurface *= 0.0;

    // find in-solid cells
    forAll(body_, cellI)
    {
        if (body_[cellI] >= 0.5)
        {
            outSurface[cellI] = 1.0;
        }
    }

    // find boundary cells
    forAll(boundaryCells_, bCell)
    {
        label cellI = boundaryCells_[bCell].first();
        outSurface[cellI] = 1.0;
    }

    //~ outSurface.write();
}

//---------------------------------------------------------------------------//
void ibInterpolation::createInnerSurface
(
    volScalarField& inSurface
)
{
    // reset field
    inSurface *= 0.0;

    // find in-solid cells
    forAll(body_, cellI)
    {
        if (body_[cellI] >= 0.5)
        {
            inSurface[cellI] = 1.0;
        }
    }

    //~ inSurface.write();
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateSurfNorm
(
)
{
    // stabilisation
    dimensionedScalar deltaN("deltaN", dimless/dimLength, SMALL);

    // calculate the surface normal based on the body gradient
    surfNorm_ = -fvc::grad(body_);
    surfNorm_ /= (mag(surfNorm_) + deltaN);
}

//---------------------------------------------------------------------------//
void ibInterpolation::correctSurfNorm
(
)
{
    // correct surface normals
    forAll(boundaryCells_, bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[bCell].first();

        // set the same surface normal to all interpolation points
        forAll(intInfoList_[bCell].intCells_, iCell)
        {
            // get the cell label
            label iCellI = intInfoList_[bCell].intCells_[iCell];

            // set the surface normal
            surfNorm_[iCellI] = surfNorm_[cellI];
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateSurfTans
(
    volVectorField& U
)
{
    // stabilisation
    dimensionedScalar deltaN("deltaN", dimLength/dimTime, SMALL);

    // reset the surface tangential field
    surfTan_ = U/(mag(U) + deltaN);

    // create interpolator
    autoPtr<interpolation<vector>> interpU = interpolation<vector>::New(HFDIBInterpDict_, U);

    // correct the tangents for boundary cells and their interpolation cells
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // get velocity in the first interpolation point
        vector UP1 = interpU->interpolate(intInfoList_[bCell].intPoints_[1], intInfoList_[bCell].intCells_[0]);

        // compute the tangential direction
        surfTan_[cellI] = UP1 - (UP1 & surfNorm_[cellI])*surfNorm_[cellI];

        // set the same surface tangent to all interpolation points
        forAll(intInfoList_[bCell].intCells_, iCell)
        {
            // get the cell label
            label iCellI = intInfoList_[bCell].intCells_[iCell];

            // set the surface tangent
            surfTan_[iCellI] = surfTan_[cellI];
        }
    }

    // unify the field
    surfTan_ /= (mag(surfTan_) + SMALL);
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateDistToBoundary
(
)
{
    forAll(boundaryCells_, bCell)
    {
        label outCellI = boundaryCells_[bCell].first();
        label inCellI = boundaryCells_[bCell].second();

        Tuple2<scalar,scalar> toSave;

        scalar ds;
        point surfPoint;

        if (body_[outCellI] >= thrSurf_)
        {
            ds = Foam::atanh(1-2*body_[outCellI])*Foam::pow(VAve_,0.333)/intSpan_; // y > 1 for lambda < 0.5
            toSave.first() = ds;

            surfPoint = mesh_.C()[outCellI];
            surfPoint -= surfNorm_[outCellI]*ds;
            ds = -1*surfNorm_[outCellI] & (mesh_.C()[inCellI] - surfPoint);
            toSave.second() = ds;
        }

        else if (body_[inCellI] < 1.0)
        {
            ds = -1*Foam::atanh(1-2*body_[inCellI])*Foam::pow(VAve_,0.333)/intSpan_; // y > 1 for lambda < 0.5
            toSave.second() = ds;

            surfPoint = mesh_.C()[inCellI];
            surfPoint += surfNorm_[inCellI]*ds;
            ds = surfNorm_[inCellI] & (mesh_.C()[outCellI] - surfPoint);
            toSave.first() = ds;

        }

        else
        {
            // find the shared faces
            label faceI;

            forAll(mesh_.cells()[inCellI], fI)
            {
                faceI = mesh_.cells()[inCellI][fI];

                if (mesh_.faceOwner()[faceI] == outCellI or mesh_.faceNeighbour()[faceI] == outCellI)
                {
                    break;
                }
            }

            // compute ds as a distance from the cell centers to the center of the shared face
            ds = mag(mesh_.C()[inCellI] - mesh_.Cf()[faceI]);
            toSave.second() = ds;

            ds = mag(mesh_.C()[outCellI] - mesh_.Cf()[faceI]);
            toSave.first() = ds;
        }

        boundaryDists_[bCell] = toSave;
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::saveInterpolationInfo
(
    word outDir,
    word fileName
)
{
    autoPtr<OFstream> outFilePtr;

    outFilePtr.reset(new OFstream(outDir/fileName));
    outFilePtr() << "cellI,cellCenter,order,intPoints,intCells" << endl;

    forAll(boundaryCells_, bCell)
    {
        outFilePtr() << boundaryCells_[bCell].first() << ","
            << mesh_.C()[boundaryCells_[bCell].first()] << ","
            << intInfoList_[bCell].order_ << ","
            << intInfoList_[bCell].intPoints_ << ","
            << intInfoList_[bCell].intCells_ << endl;
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

// ************************************************************************* //
