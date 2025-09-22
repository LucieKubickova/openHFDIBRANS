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
iProci_
(
    IOobject
    (
        "iProci",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("zero", dimless, -1.0)
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
    // prepare cell centers for acces to neighbors across processor boundary
    // Note (LK): may be redundant
    volVectorField cellCenters
    (
        IOobject
        (
            "cellCenters",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        //~ mesh_.C()
        dimensionedVector("zero", dimLength, vector::zero)
    );

    forAll(mesh_.C(), cellI)
    {
        vector center = mesh_.C()[cellI];
        cellCenters[cellI] = center;
    }

    cellCenters.correctBoundaryConditions();

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
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;
        label iFace = boundaryCells_[Pstream::myProcNo()][bCell].iFace_;

        // find surf point
        label surfCell;
        point surfPoint;
        vector surfNormToSend;
        scalar sigma;

        // separate for different kinds of boundary cells
        if (body_[outCellI] < thrSurf_)
        {
            if (Pstream::myProcNo() != iProc)
            {
                forAll(mesh_.boundaryMesh(), patchI)
                {
                    if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI]))
                    {
                        const processorPolyPatch& procPatch
                            = refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

                        // get the neighboring processor id
                        label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                            ? procPatch.neighbProcNo() : procPatch.myProcNo();

                        if (sProc == iProc)
                        {
                            // access neighbor processor patch
                            const scalarField& bodyN = body_.boundaryField()[patchI].patchNeighbourField();
                            const vectorField& surfNormN = surfNorm_.boundaryField()[patchI].patchNeighbourField();
                            const vectorField& cellCenterN = cellCenters.boundaryField()[patchI].patchNeighbourField();

                            // calc data
                            surfPoint = cellCenterN[iFace];
                            sigma = boundaryCells_[Pstream::myProcNo()][bCell].sigma_;
                            surfNormToSend = surfNormN[iFace];
                            surfPoint += surfNormN[iFace]*sigma;
                        }
                    }
                }
            }

            else
            {
                surfPoint = mesh_.C()[inCellI];;
                sigma = boundaryCells_[Pstream::myProcNo()][bCell].sigma_;
                surfNormToSend = surfNorm_[inCellI];
                surfPoint += surfNormToSend*sigma;
            }
        }
        else
        {
            surfPoint = mesh_.C()[outCellI];
            sigma = boundaryCells_[Pstream::myProcNo()][bCell].sigma_;
            surfNormToSend = surfNorm_[outCellI];
            surfPoint += surfNormToSend*sigma;
        }

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
        boundaryCells_[Pstream::myProcNo()][bCell].fProc_ = intPoints[1].iProc_;
    }

    // prepare lists for surface cells
    // Note (LK): parallelization of this was not fixed nor checked
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
void ibInterpolation::findBoundaryCells
(
)
{
    // preparation
    boundaryCell bCellToAdd;
    isBoundaryCell_ = -1;

    // loop over cells
    forAll(mesh_.cellCells(), cellI)
    {
        
        bool isBoundary(false);
        label iCell(-1);
        label iFace(-1);
        label iProc(Pstream::myProcNo());

        iProci_[cellI] = iProc; // Note (LK): save for testing

        if (body_[cellI] < 0.5 && body_[cellI] >= thrSurf_)
        {
            findNeighborInBody(cellI, thrSurf_, iCell, iFace, iProc, isBoundary);
            isBoundary = true; // Note (LK): should be always true
        }

        else if (body_[cellI] < thrSurf_)
        {
            findNeighborInBody(cellI, 0.5, iCell, iFace, iProc, isBoundary);
        }

        // check whether the cell is adjecent to a regular wall
        if (isBoundary and excludeWalls_)
        {
            isBoundary = !(isWallCell(cellI));
        }

        // add the cell
        if (isBoundary)
        {
            bCellToAdd.bCell_ = cellI;
            bCellToAdd.iCell_ = iCell;
            bCellToAdd.iProc_ = iProc;
            bCellToAdd.iFace_ = iFace;

            // save
            boundaryCells_[Pstream::myProcNo()].append(bCellToAdd);
        }
    }

    // Note (LK): save for testing
    iProci_.write();

    // sync with other processors
    List<DynamicList<label>> iFacesToSync(Pstream::nProcs());
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        label iFace = boundaryCells_[Pstream::myProcNo()][bCell].iFace_;
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;

        if (iProc != Pstream::myProcNo())
        {
            iFacesToSync[iProc].append(iFace);
        }
    }

    PstreamBuffers pBufsIFaces(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIFaces(proci, pBufsIFaces);
            sendIFaces << iFacesToSync[proci];
        }
    }

    pBufsIFaces.finishedSends();
    
    // recieve
    List<DynamicList<label>> iCellsToRetr(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIFaces(proci, pBufsIFaces);
            DynamicList<label> recIFaces (recvIFaces);

            // find cells for inner faces
            forAll(recIFaces, rFace)
            {
                // get the cell label
                label faceI = recIFaces[rFace]; // local face labels

                // find the respective cell 
                forAll(mesh_.boundaryMesh(), patchI)
                {
                    if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI]))
                    {
                        const processorPolyPatch& procPatch
                            = refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

                        // get the neighboring processor id
                        label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                            ? procPatch.neighbProcNo() : procPatch.myProcNo();

                        if (sProc == proci)
                        {
                            // get the cell label
                            label rCellI = mesh_.boundaryMesh()[patchI].faceCells()[faceI];
                            iCellsToRetr[proci].append(rCellI);
                        }
                    }
                }
            }
        }
    }

    // return
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendICells(proci, pBufsIFaces);
            sendICells << iCellsToRetr[proci];
        }
    }
    
    pBufsIFaces.finishedSends();

    List<DynamicList<label>> iCellsCmpl(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvICells(proci, pBufsIFaces);
            DynamicList<label> recICells (recvICells);
            iCellsCmpl[proci] = recICells;
        }
    }

    pBufsIFaces.clear();

    // prepare counter
    List<label> counter(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        counter[proci] = 0;
    }

    // complete boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;
        if (iProc != Pstream::myProcNo())
        {   
            // get the returned cell label
            label rCell = iCellsCmpl[iProc][counter[iProc]];
            ++counter[iProc];

            // save
            boundaryCells_[Pstream::myProcNo()][bCell].iCell_ = rCell;
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
void ibInterpolation::findNeighborInBody
(
    label& cellI,
    scalar threshold,
    label& iCell,
    label& iFace,
    label& iProc,
    bool& isBoundary
)
{
    // search
    if (boundarySearch_ == "vertex")
    {
        if (Pstream::nProcs() == 1) // Note (LK): not really working single core version, gonna be removed
        {
            forAll(mesh_.cellPoints()[cellI], pID) // vertex neighbours
            {
                label pointI = mesh_.cellPoints()[cellI][pID];
    
                forAll(mesh_.pointCells()[pointI], cI)
                {
                    if (body_[mesh_.pointCells()[pointI][cI]] >= threshold)
                    {
                        isBoundary = true;
                        iCell = mesh_.pointCells()[pointI][cI];
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
        // get labels
        label faceI = getFaceInDir(cellI, -1*surfNorm_[cellI]);
        label nI(-1);
    
        // check for non-internal cells
        if (!mesh_.isInternalFace(faceI))
        {
            // get the patch the face belongs to
            label facePatchI(mesh_.boundaryMesh().whichPatch(faceI));
            const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchI];
    
            // check if it is a processor boundary
            if (cPatch.type() == "processor")
            {
                // get the processor patch
                const processorPolyPatch& procPatch
                    = refCast<const processorPolyPatch>(cPatch);
    
                // get the neighboring processor id
                label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                    ? procPatch.neighbProcNo() : procPatch.myProcNo();
    
                // access neighbor processor patch
                const scalarField& bodyN = body_.boundaryField()[facePatchI].patchNeighbourField();
    
                // get local face value
                label localI = cPatch.whichFace(faceI);
    
                // acces neighbor value
                if (bodyN[localI] >= threshold)
                {
                    isBoundary = true;
                    iFace = localI;
                    iProc = sProc;
                }
            }
        }
    
        else
        {
            // get the neighbor cell label
            label owner(mesh_.owner()[faceI]);
            label neighbor(mesh_.neighbour()[faceI]);
            nI = (cellI == owner) ? neighbor : owner;
    
            // check lambda field
            if (body_[nI] >= threshold)
            {
                isBoundary = true;
                iCell = nI;
                iProc = Pstream::myProcNo();
            }
        }
    }
    
    else
    {
        FatalError << "Boundary cell search " << boundarySearch_ << " not implemented" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
label ibInterpolation::getFaceInDir
(
    label cellI,
    vector dir
)
{
    // prepare data
    label faceToReturn = -1;
    const labelList& cellFaces(mesh_.cells()[cellI]);

    // auxiliar scalar
    scalar dotProd(-GREAT);

    // loop over cell faces
    forAll (cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        vector outNorm = (mesh_.faceOwner()[fI] == cellI)
            ? mesh_.Sf()[fI] : (-1*mesh_.Sf()[fI]);

        scalar auxDotProd(outNorm & dir);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            faceToReturn = fI;
        }
    }

    return faceToReturn;
}

//---------------------------------------------------------------------------//
bool ibInterpolation::isWallCell
(
    label& cellI
)
{
    bool isWallCell(true);

    // get wall patches
    DynamicList<label> wPatchIs;
    forAll(mesh_.boundary(), pI)
    {
        if (mesh_.boundary()[pI].type() == "wall")
        {
            wPatchIs.append(pI);
        }
    }

    // loop over cell faces
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
                isWallCell = false;
            }
        }
    }

    return isWallCell;
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
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;

        // set the surface value
        surface[outCellI] = 0.0;

        // set inner cells
        List<DynamicList<label>> iCellsToSync(Pstream::nProcs());
        if (Pstream::myProcNo() == iProc)
        {
            surface[inCellI] = boundaryValue;
        }

        else
        {
            iCellsToSync[iProc].append(inCellI);
        }

        // sync with other processors
        PstreamBuffers pBufsICells(Pstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if(proci != Pstream::myProcNo())
            {
                UOPstream sendICells(proci, pBufsICells);
                sendICells << iCellsToSync[proci];
            }
        }

        pBufsICells.finishedSends();
        
        // recieve
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci != Pstream::myProcNo())
            {
                UIPstream recvICells(proci, pBufsICells);
                DynamicList<label> recICells (recvICells);

                // find cells for inner faces
                forAll(recICells, rCell)
                {
                    // get the cell label
                    label cellI = recICells[rCell];

                    // assign
                    surface[cellI] = boundaryValue;
                }
            }
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
    // prepare cell centers for acces to neighbors across processor boundary
    // Note (LK): may be redundant
    volVectorField cellCenters
    (
        IOobject
        (
            "cellCenters",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        //~ mesh_.C()
        dimensionedVector("zero", dimLength, vector::zero)
    );

    forAll(mesh_.C(), cellI)
    {
        vector center = mesh_.C()[cellI];
        cellCenters[cellI] = center;
    }

    cellCenters.correctBoundaryConditions();

    // loop over boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the outer and innter cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;
        label iFace = boundaryCells_[Pstream::myProcNo()][bCell].iFace_;

        // prepare
        scalar sigma(0.0);
        scalar yOrtho(0.0);
        scalar yEff(0.0);
        point surfPoint;
        scalar l(0.0);

        // if outer cell is intersected
        if (body_[outCellI] >= thrSurf_)
        {
            l = getCellSize(outCellI);
            sigma = -1*Foam::atanh(1-2*body_[outCellI])*l/intSpan_; // y < 1 for lambda < 0.5
            yOrtho = -1*sigma; // standard approach
            yEff = 0.5*(yOrtho + l*0.5);
        }

        // if inner cell is intersected
        else if (Pstream::myProcNo() == iProc)
        {
            if (body_[inCellI] < 1.0)
            {
                l = getCellSize(inCellI);
                surfPoint = mesh_.C()[inCellI];
                sigma = -1*Foam::atanh(1-2*body_[inCellI])*l/intSpan_; // y > 1 for lambda > 0.5
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
        }

        else // inner cell on a different processor
        {
            l = getCellSize(outCellI); // Note (LK): should be the size of the inner cell, needs fixing
            forAll(mesh_.boundaryMesh(), patchI)
            {
                const polyPatch& cPatch = mesh_.boundaryMesh()[patchI];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch
                        = refCast<const processorPolyPatch>(cPatch);

                    // get the neighboring processor id
                    label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                        ? procPatch.neighbProcNo() : procPatch.myProcNo();

                    if (sProc == iProc)
                    {
                        // access neighbor processor patch
                        const scalarField& bodyN = body_.boundaryField()[patchI].patchNeighbourField();
                        const vectorField& surfNormN = surfNorm_.boundaryField()[patchI].patchNeighbourField();
                        const vectorField& cellCenterN = cellCenters.boundaryField()[patchI].patchNeighbourField();

                        if (bodyN[iFace] < 1.0)
                        {
                            // calc data
                            surfPoint = cellCenterN[iFace];
                            sigma = -1*Foam::atanh(1-2*bodyN[iFace])*l/intSpan_; // y > 1 for lambda > 0.5
                            surfPoint += surfNormN[iFace]*sigma;
                            yOrtho = surfNormN[iFace] & (mesh_.C()[outCellI] - surfPoint); // standard approach

                            // effective distance
                            yEff = 0.5*(yOrtho + l*0.5);
                        }

                        // not intersected outer cells
                        else
                        {
                            if (boundarySearch_ == "face")
                            {
                                // get the processor face label
                                label faceI = cPatch.start() + iFace;

                                // get the face center
                                vector center = mesh_.Cf()[faceI];

                                // compute ds as a distance from the cell center to the averaged vertex
                                sigma = mag(cellCenterN[iFace] - center);
                                yOrtho = mag(mesh_.C()[outCellI] - center);
                                yEff = yOrtho;
                            }

                            else
                            {
                                // Note (LK): should be fixed later
                                FatalError << "Boundary cell search " << boundarySearch_ << " not implemented in parallel" << exit(FatalError);
                            }
                        }

                        break;
                    }
                }
            }
        }

        // save data
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
            //~ sigmai_[inCellI] = sigma; // Note (LK): shoud be fixed, does not work for inner cells on other processors
            sigmai_[outCellI] = sigma; // Note (LK): quick fix for checking of data, should be removed, I guess
        }
        yOrthoi_[outCellI] = yOrtho;
        yEffi_[outCellI] = yEff;
    }

    // post processing (now only possible averaging)
    postProcessYOrtho();
}

//---------------------------------------------------------------------------//
scalar ibInterpolation::getCellSize
(
    label cellI
)
{
    scalar cellSize(0.0);

    if (averageV_)
    {
        cellSize = Foam::pow(VAve_, 0.333);
    }
    else if (readL_)
    {
        cellSize = valueL_;
    }
    else
    {
        cellSize = Foam::pow(mesh_.V()[cellI], 0.333);
    }

    return cellSize;
}

//---------------------------------------------------------------------------//
void ibInterpolation::postProcessYOrtho
(
)
{
    // Note (LK): not prepared for parallel run, not usually used at all
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
        scalar l = getCellSize(cellI);

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
    DynamicList<label> saveOutCells;
    DynamicList<label> saveInCells;
    DynamicList<label> saveFreeCells;

    List<DynamicList<label>> iCellsToSync(Pstream::nProcs());
    List<DynamicList<label>> fCellsToSync(Pstream::nProcs());

    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        saveOutCells.append(boundaryCells_[Pstream::myProcNo()][bCell].bCell_);

        // check inner cells
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;
        if (Pstream::myProcNo() != iProc)
        {
            iCellsToSync[iProc].append(boundaryCells_[Pstream::myProcNo()][bCell].iCell_);
        }

        else
        {
            saveInCells.append(boundaryCells_[Pstream::myProcNo()][bCell].iCell_);
        }

        // check free stream cells
        label fProc = boundaryCells_[Pstream::myProcNo()][bCell].fProc_;
        if (Pstream::myProcNo() != fProc)
        {
            fCellsToSync[fProc].append(boundaryCells_[Pstream::myProcNo()][bCell].fCell_);
        }

        else
        {
            saveFreeCells.append(boundaryCells_[Pstream::myProcNo()][bCell].fCell_);
        }
    }

    // sync with other processors
    PstreamBuffers pBufsICells(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsFCells(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendICells(proci, pBufsICells);
            UOPstream sendFCells(proci, pBufsFCells);
            sendICells << iCellsToSync[proci];
            sendFCells << fCellsToSync[proci];
        }
    }

    pBufsICells.finishedSends();
    pBufsFCells.finishedSends();

    // recieve
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvICells(proci, pBufsICells);
            UIPstream recvFCells(proci, pBufsFCells);
            DynamicList<label> recICells (recvICells);
            DynamicList<label> recFCells (recvFCells);

            // add all inner cells 
            forAll(recICells, rCell)
            {
                saveInCells.append(recICells[rCell]);
            }

            // add all free cells
            forAll(recFCells, rCell)
            {
                saveFreeCells.append(recFCells[rCell]);
            }
        }
    }

    pBufsICells.clear();
    pBufsFCells.clear();

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

        // compute scalar product
        scalar prodUNorm = U[cellI] & surfNorm_[cellI];

        // cut U if it aims in opposite to surfNorm
        if (prodUNorm < 0.0)
        {
            U[cellI] -= (U[cellI] & surfNorm_[cellI])*surfNorm_[cellI];
        }
    }
}

//---------------------------------------------------------------------------//
void ibInterpolation::cutPhiInBoundaryCells
(
    surfaceScalarField& phi
)
{
    // Note (LK): not prepared for parallel, but not used now
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
