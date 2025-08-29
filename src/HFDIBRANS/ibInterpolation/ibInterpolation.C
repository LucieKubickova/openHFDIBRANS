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
    List<DynamicList<boundaryCell>>& boundaryCells,
    List<DynamicList<surfaceCell>>& surfaceCells,
    List<DynamicList<label>>& internalCells,
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
sigmai_
(
    IOobject
    (
        "sigmai",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("zero", dimless, -1.0)
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
yEffi_
(
    IOobject
    (
        "yEffi",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("zero", dimless, -1.0)
),
boundaryCells_(boundaryCells),
surfaceCells_(surfaceCells),
internalCells_(internalCells),
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
    boundarySearch_ = HFDIBDEMDict_.lookupOrDefault<word>("boundarySearch", "face");
    excludeWalls_ = HFDIBDEMDict_.lookupOrDefault<bool>("excludeWalls", false);
    readSurfNorm_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSurfaceNormal", false);
    aveYOrtho_ = HFDIBDEMDict_.lookupOrDefault<bool>("averageYOrtho", false);
    totalYOrthoAve_ = HFDIBDEMDict_.lookupOrDefault<bool>("totalYOrthoAverage", false);
    intSpan_ = readScalar(HFDIBDEMDict_.lookup("interfaceSpan"));
    thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
    aveCoeff_ = HFDIBDEMDict_.lookupOrDefault<scalar>("averagingCoeff", 1.0);
    nAveYOrtho_ = HFDIBDEMDict_.lookupOrDefault<label>("nAveragingYOrtho", 1.0);
    averageV_ = HFDIBDEMDict_.lookupOrDefault<bool>("averageVolume", false);
    readL_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSize", false);
    valueL_ = HFDIBDEMDict_.lookupOrDefault<scalar>("sizeValue", 0.0); // EXPERIMENTAL

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
    // prepare lists
    List<label> bCells(boundaryCells_[Pstream::myProcNo()].size());
    List<point> bPoints(boundaryCells_[Pstream::myProcNo()].size());
    List<vector> bNormals(boundaryCells_[Pstream::myProcNo()].size());

    // prepare data
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get origin cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

        // find surf point
        label surfCell;
        point surfPoint;
        scalar sigma;

        // separate for different kinds of boundary cells
        if (body_[outCellI] < thrSurf_)
        {
            surfCell = inCellI;
        }
        else
        {
            surfCell = outCellI;
        }
        surfPoint = mesh_.C()[surfCell];
        sigma = boundaryCells_[Pstream::myProcNo()][bCell].sigma_;
        vector surfNormToSend = surfNorm_[surfCell];
        surfPoint += surfNormToSend*sigma;

        // assign
        bCells[bCell] = outCellI;
        bPoints[bCell] = surfPoint;
        bNormals[bCell] = surfNormToSend;
    }

    // initialize interpolation class
    lineIntInfoBoundary_.set(new lineIntInfo(mesh_, bCells, bPoints, bNormals));

    // find interpolation points
    lineIntInfoBoundary_->setIntpInfo();

    // assign back to known structures
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        List<intPoint>& intPoints = lineIntInfoBoundary_->getIntPoints()[bCell];
        boundaryCells_[Pstream::myProcNo()][bCell].fCell_ = intPoints[1].iCell_;
    }

    // prepare lists for surface cells
    List<label> sCells(surfaceCells_[Pstream::myProcNo()].size());
    List<point> sPoints(surfaceCells_[Pstream::myProcNo()].size());
    List<vector> sNormals(surfaceCells_[Pstream::myProcNo()].size());

    // prepare data
    forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    {
        // get origin cell label
        label cellI = surfaceCells_[Pstream::myProcNo()][sCell].sCell_;

        // find surf point
        point surfPoint = mesh_.C()[cellI];
        scalar sigma = surfaceCells_[Pstream::myProcNo()][sCell].sigma_;
        vector surfNormToSend = surfNorm_[cellI];
        surfPoint += surfNormToSend*sigma;

        // assign
        sCells[sCell] = cellI;
        sPoints[sCell] = surfPoint;
        sNormals[sCell] = surfNormToSend;
    }

    // initialize interpolation class
    lineIntInfoSurface_.set(new lineIntInfo(mesh_, sCells, sPoints, sNormals));

    // find interpolation points
    lineIntInfoSurface_->setIntpInfo();
}

//---------------------------------------------------------------------------//
// Custom function to find cell containing point
// Note: this is much cheaper than standard OF functions
// Note (LK): want to remove this
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
    boundaryCell bCellToAdd;
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
                if (Pstream::nProcs() == 1) // Note (LK): working single core version, gonna be removed
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

                else
                {
                    // Note (LK): should be fixed later
                    FatalError << "Boundary cell search " << boundarySearch_ << " not implemented in parallel" << exit(FatalError);
                }
            }

            else if (boundarySearch_ == "face")
            {
                if (Pstream::nProcs() == 1) // Note (LK): working single core version, gonna be removed
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

                else // Note (LK): parallel version of boundary search face
                {
                    FatalError << "Boundary cell search " << boundarySearch_ << " not implemented in parallel" << exit(FatalError);
                    //~ forAll(mesh_.cells()[cellI], fI)
                    //~ {
                        //~ // get labels
                        //~ label faceI = mesh_.cells()[cellI][fI];
                        //~ label nI(-1);

                        //~ // check for non-internal cells
                        //~ if (!mesh_.isInternalFace(faceI))
                        //~ {
                            //~ // get the patch the face belongs to
                            //~ label facePatchI(mesh_.boundaryMesh().whichPatch(faceI));
                            //~ const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchI];

                            //~ // check if it is a processor boundary
                            //~ if (cPatch.type() == "processor")
                            //~ {
                                //~ // get the processor patch
                                //~ const processorPolyPatch& procPatch
                                    //~ = refCast<const processorPolyPatch>(cPatch);

                                //~ // get the neighboring processor id
                                //~ label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                                    //~ ? procPatch.neighbProcNo() : procPatch.myProcNo();

                                //~ // save labels
                                //~ nI = cPatch.whichFace(faceInDir);

                                //~ // check lambda
                                //~ //~ retP.iProc_ = sProc; // Note (LK): save the processor label too

                                //~ //~ return retP;
                            //~ }
                            //~ else
                            //~ {
                                //~ //~ retP.iProc_ = -1; // Note (LK): save the processor label too
                                //~ //~ return retP;
                            //~ }
                        //~ }

                        //~ // get the neighbor cell label
                        //~ label owner(mesh_.owner()[faceI]);
                        //~ label neighbor(mesh_.neighbour()[faceI]);
                        //~ nI = (cellI == owner) ? neighbor : owner;

                        //~ // check lambda field
                        //~ if (body_[nI] >= 0.5)
                        //~ {
                            //~ toInclude = true;
                            //~ vertex = nI;
                            //~ break;
                        //~ }
                    //~ }
                }
            }

            else
            {
                FatalError << "Boundary cell search " << boundarySearch_ << " not implemented" << exit(FatalError);
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

                bCellToAdd.bCell_ = cellI;
                bCellToAdd.iCell_ = startCell.second().first();
            }

            else
            {
                bCellToAdd.bCell_ = cellI;
                bCellToAdd.iCell_ = vertex;
            }

            // fix non-existent inner cells
            if (bCellToAdd.iCell_ == -1)
            {
                bCellToAdd.iCell_ = bCellToAdd.bCell_;
            }

            // assign dummy as third (first free stream cell)
            bCellToAdd.fCell_ = -1;
            boundaryCells_[Pstream::myProcNo()].append(bCellToAdd);
        }
    }

    // prepare label field
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

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
        surfaceCell sCellToAdd;

        // check lambda
        if (body_[cellI] >= thrSurf_)
        {
            if (body_[cellI] <= (1-thrSurf_))
            {
                sCellToAdd.sCell_ = cellI;
                surfaceCells_[Pstream::myProcNo()].append(sCellToAdd);
            }

            else
            {
                internalCells_[Pstream::myProcNo()].append(cellI);
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
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

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
void ibInterpolation::setOnlyInnerSurface
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
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

        // set the surface value
        surface[outCellI] = 0.0;
        surface[inCellI] = boundaryValue;
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
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

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
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the outer and innter cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

        // prepare
        scalar sigma;
        scalar yOrtho;
        scalar yEff;
        point surfPoint;
        scalar l;

        // if outer cell is intersected
        if (body_[outCellI] >= thrSurf_)
        {
            if (averageV_)
            {
                l = Foam::pow(VAve_, 0.333);
            }
            else if (readL_)
            {
                l = valueL_;
            }
            else
            {
                l = Foam::pow(mesh_.V()[outCellI], 0.333);
            }
            sigma = -1*Foam::atanh(1-2*body_[outCellI])*l/intSpan_; // y < 1 for lambda < 0.5
            yOrtho = -1*sigma; // standard approach
            yEff = 0.5*(yOrtho + l*0.5);
        }

        // if inner cell is intersected
        else if (body_[inCellI] < 1.0)
        {
            if (averageV_)
            {
                l = Foam::pow(VAve_, 0.333);
            }
            else if (readL_)
            {
                l = valueL_;
            }
            else
            {
                l = Foam::pow(mesh_.V()[inCellI], 0.333);
            }
            sigma = -1*Foam::atanh(1-2*body_[inCellI])*l/intSpan_; // y > 1 for lambda > 0.5
            surfPoint = mesh_.C()[inCellI];
            surfPoint += surfNorm_[inCellI]*sigma;
            yOrtho = surfNorm_[inCellI] & (mesh_.C()[outCellI] - surfPoint); // standard approach
            yEff = 0.5*(yOrtho + l*0.5);
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
            yEff = yOrtho;
        }

        boundaryCells_[Pstream::myProcNo()][bCell].sigma_ = sigma;
        boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_ = yOrtho;
        boundaryCells_[Pstream::myProcNo()][bCell].yEff_ = yEff;

        // save to fields
        if (body_[outCellI] >= thrSurf_)
        {
            sigmai_[outCellI] = sigma;
        }
        else
        {
            sigmai_[inCellI] = sigma;
        }
        yOrthoi_[outCellI] = yOrtho;
        yEffi_[outCellI] = yEff;
    }

    // average orthogonal distance over all boundary cells
    if (totalYOrthoAve_)
    {
        // preparation
        scalar yAve(0.0);

        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            yAve += boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
        }

        // calculate average
        yAve /= boundaryCells_[Pstream::myProcNo()].size();

        // save to all boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_ = yAve;
        }
    }

    // average orthogonal distance over neighbors
    else if (aveYOrtho_)
    {
        for (label nCorr = 0; nCorr < nAveYOrtho_; nCorr++)
        {
            // loop over boundary cells
            forAll(boundaryCells_[Pstream::myProcNo()], bCell)
            {
                // get the cell label
                label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

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
                        yAve += dotNorm/delta*boundaryCells_[Pstream::myProcNo()][nCell].yOrtho_;
                        scalesSum += dotNorm/delta;
                        nAdded += 1;
                    }
                }

                // compute average scale
                scalar aveScale = scalesSum/nAdded;

                // add current cell
                yAve += aveCoeff_*aveScale*boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
                scalesSum += aveCoeff_*aveScale;

                // divide average y by scales sum
                yAve /= scalesSum;

                // assign
                boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_ = yAve;
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
    forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    {
        // get cell label
        label cellI = surfaceCells_[Pstream::myProcNo()][sCell].sCell_;

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
        sigma = -1*Foam::atanh(1-2*body_[cellI])*l/intSpan_;

        // assign
        surfaceCells_[Pstream::myProcNo()][sCell].sigma_ = sigma;
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
    outFilePtr() << "cellI,inCellI,freeCellI,outCellCenter,inCellCenter,freeCellCenter,surfNorm,surfNormIn,bodyOut,bodyIn,sigma,yOrtho,yEff,order,intPoints,intCells" << endl;

    // Note (LK): this needs fixing 
    // loop over boundary cells
    //~ forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    //~ {
        //~ outFilePtr() << boundaryCells_[Pstream::myProcNo()][bCell].bCell_ << ","
            //~ << boundaryCells_[Pstream::myProcNo()][bCell].iCell_ << ","
            //~ << boundaryCells_[Pstream::myProcNo()][bCell].fCell_ << ","
            //~ << mesh_.C()[boundaryCells_[Pstream::myProcNo()][bCell].bCell_] << ","
            //~ << mesh_.C()[boundaryCells_[Pstream::myProcNo()][bCell].iCell_] << ","
            //~ << mesh_.C()[boundaryCells_[Pstream::myProcNo()][bCell].fCell_] << ","
            //~ << surfNorm_[boundaryCells_[Pstream::myProcNo()][bCell].bCell_] << ","
            //~ << surfNorm_[boundaryCells_[Pstream::myProcNo()][bCell].iCell_] << ","
            //~ << body_[boundaryCells_[Pstream::myProcNo()][bCell].bCell_] << ","
            //~ << body_[boundaryCells_[Pstream::myProcNo()][bCell].iCell_] << ","
            //~ << boundaryCells_[Pstream::myProcNo()][bCell].sigma_ << ","
            //~ << boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_ << ","
            //~ << boundaryCells_[Pstream::myProcNo()][bCell].yEff_ << ","
            //~ << intInfoListBoundary_[Pstream::myProcNo()][bCell].order_ << ","
            //~ << intInfoListBoundary_[Pstream::myProcNo()][bCell].intPoints_ << ","
            //~ << intInfoListBoundary_[Pstream::myProcNo()][bCell].intCells_ << endl;
    //~ }

    // prepare file
    fileName = name + "_surface.dat";
    outFilePtr.reset(new OFstream(outDir/fileName));
    outFilePtr() << "cellI,cellCenter,surfNorm,body,sigma,order,intPoints,intCells" << endl;

    //~ // loop over surface cells
    //~ forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    //~ {
        //~ outFilePtr () << surfaceCells_[Pstream::myProcNo()][sCell].sCell_ << ","
            //~ << mesh_.C()[surfaceCells_[Pstream::myProcNo()][sCell].sCell_] << ","
            //~ << surfNorm_[surfaceCells_[Pstream::myProcNo()][sCell].sCell_] << ","
            //~ << body_[surfaceCells_[Pstream::myProcNo()][sCell].sCell_] << ","
            //~ << surfaceCells_[Pstream::myProcNo()][sCell].sigma_ << ","
            //~ << intInfoListSurface_[Pstream::myProcNo()][sCell].order_ << ","
            //~ << intInfoListSurface_[Pstream::myProcNo()][sCell].intPoints_ << ","
            //~ << intInfoListSurface_[Pstream::myProcNo()][sCell].intCells_ << endl;
    //~ }
}

//---------------------------------------------------------------------------//
void ibInterpolation::saveBoundaryCells
(
)
{
    List<label> saveOutCells(boundaryCells_[Pstream::myProcNo()].size());
    List<label> saveInCells(boundaryCells_[Pstream::myProcNo()].size());
    List<label> saveFreeCells(boundaryCells_[Pstream::myProcNo()].size());

    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        saveOutCells[bCell] = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        saveInCells[bCell] = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;
        saveFreeCells[bCell] = boundaryCells_[Pstream::myProcNo()][bCell].fCell_;
    }

    saveCellSet(saveOutCells, "outerBoundaryCells");
    saveCellSet(saveInCells, "innerBoundaryCells");
    saveCellSet(saveFreeCells, "freeBoundaryCells");
}

//---------------------------------------------------------------------------//
void ibInterpolation::saveSurfaceCells
(
)
{
    List<label> saveCells(surfaceCells_[Pstream::myProcNo()].size());

    forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    {
        saveCells[sCell] = surfaceCells_[Pstream::myProcNo()][sCell].sCell_;
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
    // prepare cell set
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

    // get name of output folder
    word outDir;
    if (Pstream::nProcs() == 1)
    {
        outDir = "./constant/polyMesh/sets";
    }
    else
    {
        outDir = "./processor" + Foam::name(Pstream::myProcNo()) + "/constant/polyMesh/sets";
    }

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
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

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
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell labels
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        //~ label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        //~ label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

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
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get the cell label
            label cellB = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

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
