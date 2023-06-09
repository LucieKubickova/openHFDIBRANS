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
    List<Tuple2<scalar,scalar>>& boundaryDists
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
boundaryCells_(boundaryCells),
boundaryDists_(boundaryDists),
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

            // check whether the cell is adjecent to a regular wall
            if (toInclude)
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
    // stabilisation
    dimensionedScalar deltaN("deltaN", dimless/dimLength, SMALL);

    // calculate the surface normal based on the body gradient
    surfNorm_ = -fvc::grad(body_);
    surfNorm_ /= (mag(surfNorm_) + deltaN);
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
            ds = mag(mesh_.C()[inCellI] - center);
            toSave.second() = ds;

            ds = mag(mesh_.C()[outCellI] - center);
            toSave.first() = ds;
        }

        boundaryDists_[bCell] = toSave;
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
