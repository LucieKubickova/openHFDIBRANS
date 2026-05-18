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
    ibMesh& ibMesh,
    const volScalarField& body,
    List<DynamicList<boundaryCell>>& boundaryCells,
    List<DynamicList<surfaceCell>>& surfaceCells,
    List<DynamicList<label>>& internalCells,
    labelField& isBoundaryCell
)
:
mesh_(mesh),
ibMesh_(ibMesh),
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
    dimensionedVector("zero", dimless, vector::zero)
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
    excludePatch_ = HFDIBDEMDict_.lookupOrDefault<word>("excludePatch", "none");
    readSurfNorm_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSurfaceNormal", false);
    aveYOrtho_ = HFDIBDEMDict_.lookupOrDefault<bool>("averageYOrtho", false);
    totalYOrthoAve_ = HFDIBDEMDict_.lookupOrDefault<bool>("totalYOrthoAverage", false);
    intSpan_ = readScalar(HFDIBDEMDict_.lookup("interfaceSpan"));
    thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
    aveCoeff_ = HFDIBDEMDict_.lookupOrDefault<scalar>("averagingCoeff", 1.0);
    nAveYOrtho_ = HFDIBDEMDict_.lookupOrDefault<label>("nAveragingYOrtho", 1.0);
    innerThres_ = HFDIBDEMDict_.lookupOrDefault<scalar>("innerThreshold", 0.5); // LK: experimental
    sdBasedLambda_ = HFDIBDEMDict_.lookupOrDefault<bool>("sdBasedLambda", true);
    correctIntPoints_ = HFDIBDEMDict_.lookupOrDefault<bool>("correctIntPoints", false);
    yFromCutEdges_ = HFDIBDEMDict_.lookupOrDefault<bool>("yFromCutEdges", false); // LK: experimental

    // read fvSchemes
    HFDIBInnerSchemes_ = fvSchemes_.subDict("HFDIBSchemes").subDict("innerSchemes");

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

        // find surface point
        point surfPoint(vector::zero);
        vector surfNormToSend(vector::zero);
        scalar sigma(0.0);

        // separate for different kinds of boundary cells
        if (sdBasedLambda_)
        {
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
                                const tmp<vectorField> tsurfNormN(surfNorm_.boundaryField()[patchI].patchNeighbourField());
                                const tmp<vectorField> tcellCenterN(cellCenters.boundaryField()[patchI].patchNeighbourField());
                                const vectorField& surfNormN = tsurfNormN();
                                const vectorField& cellCenterN = tcellCenterN();

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
                    surfPoint = mesh_.C()[inCellI];
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

            // save for later
            boundaryCells_[Pstream::myProcNo()][bCell].sPoint_ = surfPoint;
        }

        else
        {
            // Note (LK): data calculated in boundary dist calculation from stl file
            surfPoint = boundaryCells_[Pstream::myProcNo()][bCell].sPoint_;
            surfNormToSend = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;
        }

        // assign
        bCells[bCell] = outCellI;
        bPoints[bCell] = surfPoint;
        bNormals[bCell] = surfNormToSend;
    }

    // initialize interpolation class
    lineIntInfoBoundary_.set(new lineIntInfo(mesh_, ibMesh_, bCells, bPoints, bNormals, correctIntPoints_));

    // find interpolation points
    lineIntInfoBoundary_->setIntpInfo();

    // assign back to known structures
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        List<intPoint>& intPoints = lineIntInfoBoundary_->getIntPoints()[bCell];
        if (!intPoints[0].last_)
        {
            boundaryCells_[Pstream::myProcNo()][bCell].fPoint1_ = intPoints[1].iPoint_;
            boundaryCells_[Pstream::myProcNo()][bCell].fCell1_ = intPoints[1].iCell_; 
            boundaryCells_[Pstream::myProcNo()][bCell].fProc1_ = intPoints[1].iProc_;
        }

        if (!intPoints[1].last_)
        {
            boundaryCells_[Pstream::myProcNo()][bCell].fPoint2_ = intPoints[2].iPoint_;
            boundaryCells_[Pstream::myProcNo()][bCell].fCell2_ = intPoints[2].iCell_; 
            boundaryCells_[Pstream::myProcNo()][bCell].fProc2_ = intPoints[2].iProc_;
        }
    }

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
    lineIntInfoSurface_.set(new lineIntInfo(mesh_, ibMesh_, sCells, sPoints, sNormals, correctIntPoints_));

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
            findNeighborInBody(cellI, innerThres_, iCell, iFace, iProc, isBoundary);
        }

        // check whether the cell is adjecent to a regular wall
        if (isBoundary and excludeWalls_)
        {
            isBoundary = !(ibMesh_.isWallCell(cellI));
        }

        else if (isBoundary and excludePatch_ != "none")
        {
            isBoundary = !(ibMesh_.isOnPatch(cellI, excludePatch_));
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
    //~ iProci_.write();

    // sync with other processors
    List<DynamicList<label>> iFacesToSync(Pstream::nProcs());
    List<DynamicList<label>> bLabelsToRecv(Pstream::nProcs());

    // loop over boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        label iFace = boundaryCells_[Pstream::myProcNo()][bCell].iFace_;
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;

        if (iProc != Pstream::myProcNo())
        {
            iFacesToSync[iProc].append(iFace);
            bLabelsToRecv[iProc].append(bCell);
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
        if (proci != Pstream::myProcNo())
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

    // complete boundary cells
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            forAll(iCellsCmpl[proci], rCell)
            {
                label iCell = iCellsCmpl[proci][rCell];
                label bCell = bLabelsToRecv[proci][rCell];
                boundaryCells_[Pstream::myProcNo()][bCell].iCell_ = iCell;
            }
        }
    }

    // prepare label field
    isBoundaryCell_ *= 0.0;
    isBoundaryCell_ += -1;
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // assign
        isBoundaryCell_[cellI] = bCell;
    }

    // assign surface normal to boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // save surface normal
        boundaryCells_[Pstream::myProcNo()][bCell].sNorm_ = surfNorm_[cellI];
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

    else if (boundarySearch_ == "edge")
    {
        // get the best face, edge and vertex
        vector dir = -1*surfNorm_[cellI];
        label faceI = ibMesh_.getFaceInDir(cellI, dir);
        label edgeI = ibMesh_.getEdgeInDir(faceI, cellI, dir);
        //~ label vertI = ibMesh_.getVertInDir(edgeI, cellI, dir);

        label faceNI(-1);
        label edgeNI(-1);
    
        // check for non-internal cells
        if (!mesh_.isInternalFace(faceI))
        {
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
    
                //~ // access neighbor processor patch
                //~ const scalarField& bodyN = body_.boundaryField()[facePatchI].patchNeighbourField();
    
                //~ // get local face value
                //~ label localI = cPatch.whichFace(faceI);
    
                //~ // acces neighbor value
                //~ if (bodyN[localI] >= threshold)
                //~ {
                    //~ isBoundary = true;
                    //~ iFace = localI;
                    //~ iProc = sProc;
                //~ }
            //~ }
        }
    
        else
        {
            // get the face neighbor
            label owner(mesh_.owner()[faceI]);
            label neighbor(mesh_.neighbour()[faceI]);
            faceNI = (cellI == owner) ? neighbor : owner;

            // get the edge neighbor
            bool found(false);
            forAll(mesh_.edgeCells()[edgeI], cI)
            {
                // get cell label
                edgeNI = mesh_.edgeCells()[edgeI][cI];

                // skip current cell and the face neighbor
                if (edgeNI == cellI or edgeNI == faceNI)
                {
                    continue;
                }

                // check if vertex neighbor is neighbor with face neighbor O:)
                forAll(mesh_.cellCells()[faceNI], pI)
                {
                    label posNI = mesh_.cellCells()[faceNI][pI];
                    if (posNI == edgeNI)
                    {
                        found = true;
                    }
                }

                if (found)
                {
                    break;
                }
            }
    
            // check lambda field
            // LK: needs check if edgeNI found
            // LK: needs check which more in dir
            if (body_[faceNI] >= threshold)
            {
                isBoundary = true;
                iCell = faceNI;
                iProc = Pstream::myProcNo();
            }

            else if (body_[edgeNI] >= threshold)
            {
                isBoundary = true;
                iCell = edgeNI;
                iProc = Pstream::myProcNo();
            }
        }

        // Note (LK): should be fixed later
        FatalError << "Boundary cell search " << boundarySearch_ << " not finished" << exit(FatalError);
    }
    
    else if (boundarySearch_ == "face")
    {
        // get labels
        // Note (LK): surf norm should be from STL file if stated
        vector dir = -1*surfNorm_[cellI];
        label faceI = ibMesh_.getFaceInDir(cellI, dir);
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
                const tmp<scalarField> tbodyN(body_.boundaryField()[facePatchI].patchNeighbourField());
                const scalarField& bodyN = tbodyN();
    
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
void ibInterpolation::correctSurfaceByNormal
(
    volSymmTensorField& normSurface,
    volScalarField& surface,
    scalar bodyOnLimit
)
{
    // loop over cells
    forAll(mesh_.C(), cellI)
    {
        // check scalar surface
        if (surface[cellI] == 0.0)
        {
            continue;
        }

        // check body
        if (body_[cellI] >= bodyOnLimit)
        {
            normSurface[cellI].xx() = 1.0;
            normSurface[cellI].yy() = 1.0;
            normSurface[cellI].zz() = 1.0;
            continue;
        }

        // correct surface
        normSurface[cellI].xx() = mag(surface[cellI]*surfNorm_[cellI].x());
        normSurface[cellI].yy() = mag(surface[cellI]*surfNorm_[cellI].y());
        normSurface[cellI].zz() = mag(surface[cellI]*surfNorm_[cellI].z());
    }

    //~ // correct inner boundary cells 
    //~ forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    //~ {
        //~ // get the cell label
        //~ label iCell = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;
        //~ label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;

        //~ if (iProc == Pstream::myProcNo())
        //~ {
            //~ normSurface[iCell].xx() = mag(surface[iCell]*surfNorm_[iCell].x());
            //~ normSurface[iCell].yy() = mag(surface[iCell]*surfNorm_[iCell].y());
            //~ normSurface[iCell].zz() = mag(surface[iCell]*surfNorm_[iCell].z());
        //~ }

        //~ else
        //~ {
            //~ // Note (LK): should be fixed later
            //~ FatalError << "Surface correction by normal not implemented in parallel" << exit(FatalError);
        //~ }
    //~ }
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateSurfNorm
(
)
{
    if (readSurfNorm_)
    {
        // pass
    }

    // calculate body gradient
    volVectorField gradBody = fvc::grad(body_);

    // convert to dimmless
    gradBody *= dimensionedScalar("one", dimLength, 1.0);

    // calculate surface normal
    surfNorm_ = -gradBody;
    surfNorm_ /= (mag(surfNorm_) + SMALL);
}

//---------------------------------------------------------------------------//
void ibInterpolation::calculateSurfArea
(
)
{
    // Note (LK): not prepared for parallel
    if (Pstream::nProcs() != 1)
    {
        FatalError << "Surface area calculation not implemented in parallel" << exit(FatalError);
    }

    // Note (LK): check total area
    volScalarField surfArea
    (
        IOobject
        (
            "surfArea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    scalar totalArea(0.0);

    // Note (LK): new version
    // first surface cells
    forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    {
        // get cell label
        label cellI = surfaceCells_[Pstream::myProcNo()][sCell].sCell_;

        // prepare surf point and normal
        vector normal = surfaceCells_[Pstream::myProcNo()][sCell].sNorm_;
        point surfPoint = surfaceCells_[Pstream::myProcNo()][sCell].sPoint_;

        if (sdBasedLambda_)
        {
            normal = surfNorm_[cellI];
            scalar sigma = surfaceCells_[Pstream::myProcNo()][sCell].sigma_;
            surfPoint = mesh_.C()[cellI] + sigma*normal;

            // Note (LK): emergency assign
            surfaceCells_[Pstream::myProcNo()][sCell].sNorm_ = normal;
            surfaceCells_[Pstream::myProcNo()][sCell].sPoint_ = surfPoint;
        }

        // construct cut cell
        scalar sArea(0.0);
        sArea = ibMesh_.createCutCellAndSurface(cellI, normal, surfPoint);

        // Note (LK): check surface areas
        surfArea[cellI] = sArea;
        totalArea += sArea;

        // assign
        surfaceCells_[Pstream::myProcNo()][sCell].sArea_ = sArea;
    }

    // now boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell labels
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

        // get surface data
        vector surfNorm = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;
        point surfPoint = boundaryCells_[Pstream::myProcNo()][bCell].sPoint_;

        // prepare
        scalar sArea(0.0);

        // check body
        if (body_[outCellI] < 0.5 && body_[outCellI] >= thrSurf_)
        {
            // check if done in surface cells, else do cut
            if (surfArea[outCellI] > 0.0)
            {
                sArea = surfArea[outCellI];
            }

            else if (surfArea[inCellI] > 0.0)
            {
                sArea = surfArea[inCellI];
            }

            else
            {
                // create cut cell
                sArea = ibMesh_.createCutCellAndSurface(inCellI, surfNorm, surfPoint);
                surfArea[inCellI] = sArea;
                totalArea += sArea;
            }
        }

        else if (body_[outCellI] < thrSurf_ && body_[inCellI] < 1.0 - thrSurf_)
        {
            // check if done in surface cells, else do cut
            if (surfArea[inCellI] > 0.0)
            {
                sArea = surfArea[inCellI];
            }

            else
            {
                // create cut cell
                sArea = ibMesh_.createCutCellAndSurface(inCellI, surfNorm, surfPoint);
                surfArea[inCellI] = sArea;
                totalArea += sArea;
            }
        }

        else
        {
            // find the shared face
            forAll(mesh_.cells()[outCellI], fI)
            {
                // get face label
                label faceI = mesh_.cells()[outCellI][fI];
                
                // get owner and neighbor
                label owner(mesh_.owner()[faceI]);
                label neighbor(mesh_.neighbour()[faceI]);

                // check if shared
                if ((outCellI == owner and inCellI == neighbor) or (outCellI == neighbor and inCellI == owner))
                {
                    sArea = mag(mesh_.Sf()[faceI]);
                    break;
                }
            }

            surfArea[outCellI] = sArea;
            totalArea += sArea;
        }

        // assign
        boundaryCells_[Pstream::myProcNo()][bCell].sArea_ = sArea;
    }

    // Note (LK): check write
    surfArea.write();

    Info << "Total area is:" << endl;
    Info << "    S = " << totalArea << endl;
    Info << endl;
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
        // get the outer and inner cell label
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
            l = ibMesh_.getCellSize(outCellI);

            if (sdBasedLambda_)
            {
                sigma = -1*Foam::atanh(1-2*body_[outCellI])*l/intSpan_; // y < 0 for lambda < 0.5
            }

            else if (yFromCutEdges_)
            {
                vector surfNorm(vector::zero);
                ibMesh_.createCutCellAndCenter(outCellI, surfNorm, surfPoint);

                sigma = -1*mag(surfPoint - mesh_.C()[outCellI]);
                surfNorm /= mag(surfNorm);

                // Note (LK): check for inward pointing surface normal
                scalar dotProd = surfNorm & surfNorm_[outCellI];
                if (dotProd < 0.0)
                {
                    surfNorm *= -1;
                }

                surfNorm_[outCellI] = surfNorm/mag(surfNorm);

                // Note (LK): save to boundary cell, has to be done somewhere centraly
                boundaryCells_[Pstream::myProcNo()][bCell].sNorm_ = surfNorm_[outCellI];
                boundaryCells_[Pstream::myProcNo()][bCell].sPoint_ = surfPoint;
            }

            else
            {
                scalar intDist = Foam::pow(mesh_.V()[outCellI],0.333);
                vector surfPoint = vector::zero;
                vector surfNorm = vector::zero;
                
                ibMesh_.getClosestPointAndNormal(
                    mesh_.C()[outCellI],
                    intDist*2*vector::one,
                    surfPoint,
                    surfNorm 
                );

                sigma = -1*mag(surfPoint - mesh_.C()[outCellI]);
                surfNorm_[outCellI] = surfNorm/mag(surfNorm);

                // Note (LK): save to boundary cell, has to be done somewhere centraly
                boundaryCells_[Pstream::myProcNo()][bCell].sNorm_ = surfNorm_[outCellI];
                boundaryCells_[Pstream::myProcNo()][bCell].sPoint_ = surfPoint;
            }

            // compute distances
            yOrtho = -1*sigma; // standard approach
            yEff = 0.5*(yOrtho + l*0.5);
        }

        // if inner cell is intersected
        else if (Pstream::myProcNo() == iProc)
        {
            if (body_[inCellI] < 1.0 - thrSurf_ and sdBasedLambda_)
            {
                l = ibMesh_.getCellSize(inCellI);
                surfPoint = mesh_.C()[inCellI];
                sigma = -1*Foam::atanh(1-2*body_[inCellI])*l/intSpan_; // y > 0 for lambda > 0.5
                surfPoint += surfNorm_[inCellI]*sigma;
                yOrtho = surfNorm_[inCellI] & (mesh_.C()[outCellI] - surfPoint); // standard approach
                yEff = 0.5*(yOrtho + l*0.5);
            }

            else if (body_[inCellI] < 1.0 - thrSurf_ and yFromCutEdges_)
            {
                vector surfNorm(vector::zero);
                ibMesh_.createCutCellAndCenter(inCellI, surfNorm, surfPoint);
                surfNorm /= mag(surfNorm);

                // Note (LK): check for inward pointing surface normal
                scalar dotProd = surfNorm & surfNorm_[outCellI];
                if (dotProd < 0.0)
                {
                    surfNorm *= -1;
                }

                sigma = mag(surfPoint - mesh_.C()[inCellI]);
                surfNorm_[inCellI] = surfNorm/mag(surfNorm);
                yOrtho = surfNorm_[inCellI] & (mesh_.C()[outCellI] - surfPoint); // standard approach
                l = ibMesh_.getCellSize(outCellI);
                yEff = 0.5*(yOrtho + l*0.5);

                // Note (LK): save to boundary cell, has to be done somewhere centraly
                boundaryCells_[Pstream::myProcNo()][bCell].sNorm_ = surfNorm_[inCellI];
                boundaryCells_[Pstream::myProcNo()][bCell].sPoint_ = surfPoint;
            }

            // not intersected outer cells
            else if (sdBasedLambda_)
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

            else
            {
                scalar intDist = Foam::pow(mesh_.V()[outCellI],0.333);
                vector surfPoint = vector::zero;
                vector surfNorm = vector::zero;

                ibMesh_.getClosestPointAndNormal(
                    mesh_.C()[outCellI],
                    intDist*2*vector::one,
                    surfPoint,
                    surfNorm 
                );

                sigma = mag(surfPoint - mesh_.C()[inCellI]); // LK: not really true
                yOrtho = mag(surfPoint - mesh_.C()[outCellI]);
                surfNorm_[outCellI] = surfNorm/mag(surfNorm);

                // Note (LK): save to boundary cell, has to be done somewhere centraly
                boundaryCells_[Pstream::myProcNo()][bCell].sNorm_ = surfNorm_[outCellI];
                boundaryCells_[Pstream::myProcNo()][bCell].sPoint_ = surfPoint;

                l = ibMesh_.getCellSize(outCellI);
                yEff = 0.5*(yOrtho + l*0.5);
            }
        }

        else // inner cell on a different processor
        {
            if (sdBasedLambda_)
            {
                l = ibMesh_.getCellSize(outCellI); // Note (LK): should be the size of the inner cell, needs fixing
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
                            const tmp<scalarField> tbodyN(body_.boundaryField()[patchI].patchNeighbourField());
                            const tmp<vectorField> tsurfNormN(surfNorm_.boundaryField()[patchI].patchNeighbourField());
                            const tmp<vectorField> tcellCenterN(cellCenters.boundaryField()[patchI].patchNeighbourField());

                            const scalarField& bodyN = tbodyN();
                            const vectorField& surfNormN = tsurfNormN();
                            const vectorField& cellCenterN = tcellCenterN();

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

            else
            {
                scalar intDist = Foam::pow(mesh_.V()[outCellI],0.333);
                vector surfPoint = vector::zero;
                vector surfNorm = vector::zero;

                ibMesh_.getClosestPointAndNormal(
                    mesh_.C()[outCellI],
                    intDist*2*vector::one,
                    surfPoint,
                    surfNorm 
                );

                sigma = mag(surfPoint - mesh_.C()[inCellI]); // LK: not really true
                yOrtho = mag(surfPoint - mesh_.C()[outCellI]);
                surfNorm_[outCellI] = surfNorm/mag(surfNorm);

                // Note (LK): save to boundary cell, has to be done somewhere centraly
                boundaryCells_[Pstream::myProcNo()][bCell].sNorm_ = surfNorm_[outCellI];
                boundaryCells_[Pstream::myProcNo()][bCell].sPoint_ = surfPoint;

                l = ibMesh_.getCellSize(outCellI);
                yEff = 0.5*(yOrtho + l*0.5);
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

    // Note (LK): write check, remove
    surfNorm_.write();
    sigmai_.write();
    yOrthoi_.write();
    yEffi_.write();

    // post processing (now only possible averaging)
    postProcessYOrtho();
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
        scalar l = ibMesh_.getCellSize(cellI);

        // calculate signed distance
        if (sdBasedLambda_)
        {
            sigma = -1*Foam::atanh(1-2*body_[cellI])*l/intSpan_;
        }

        else
        {
            scalar intDist = Foam::pow(mesh_.V()[cellI],0.333);
            vector surfPoint = vector::zero;
            vector surfNorm = vector::zero;

            ibMesh_.getClosestPointAndNormal(
                mesh_.C()[cellI],
                intDist*2*vector::one,
                surfPoint,
                surfNorm
            );

            sigma = mag(surfPoint - mesh_.C()[cellI]);
            if (body_[cellI] < 0.5)
            {
                sigma *= -1;
            }
            surfNorm_[cellI] = surfNorm/mag(surfNorm);

            // Note (LK): data calculated in boundary dist calculation from stl file
            surfaceCells_[Pstream::myProcNo()][sCell].sNorm_ = surfNorm_[cellI];
            surfaceCells_[Pstream::myProcNo()][sCell].sPoint_ = surfPoint;
        }

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
    outFilePtr() << "cellI" << ","
        << "surfNorm_x" << ","
        << "surfNorm_y" << ","
        << "surfNorm_z" << ","
        << "sigma" << ","
        << "yOrtho" << ","
        << "intCell0" << ","
        << "intCell1" << ","
        << "intCell2" << ","
        << "intPoint0_x" << ","
        << "intPoint0_y" << ","
        << "intPoint0_z" << ","
        << "intPoint1_x" << ","
        << "intPoint1_y" << ","
        << "intPoint1_z" << ","
        << "intPoint2_x" << ","
        << "intPoint2_y" << ","
        << "intPoint2_z" << endl;

    // loop over boundary cells
    if (Pstream::nProcs() == 1)
    {
        List<List<intPoint>>& intInfos = lineIntInfoBoundary_->getIntPoints();

        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            boundaryCell cellToSave = boundaryCells_[Pstream::myProcNo()][bCell];

            outFilePtr() << cellToSave.bCell_ << ","
                << surfNorm_[cellToSave.bCell_].x() << ","
                << surfNorm_[cellToSave.bCell_].y() << ","
                << surfNorm_[cellToSave.bCell_].z() << ","
                << cellToSave.sigma_ << ","
                << cellToSave.yOrtho_ << ","
                << intInfos[bCell][0].iCell_ << ","
                << intInfos[bCell][1].iCell_ << ","
                << intInfos[bCell][2].iCell_ << ","
                << intInfos[bCell][0].iPoint_.x() << ","
                << intInfos[bCell][0].iPoint_.y() << ","
                << intInfos[bCell][0].iPoint_.z() << ","
                << intInfos[bCell][1].iPoint_.x() << ","
                << intInfos[bCell][1].iPoint_.y() << ","
                << intInfos[bCell][1].iPoint_.z() << ","
                << intInfos[bCell][2].iPoint_.x() << ","
                << intInfos[bCell][2].iPoint_.y() << ","
                << intInfos[bCell][2].iPoint_.z() << endl;
        }
    }


    // Note (LK): OLD SAVE SYSTEM
    //~ outFilePtr() << "cellI,inCellI,freeCellI1,freeCellI2,outCellCenter,inCellCenter,freeCellCenter1,freeCellCenter2,surfNorm,surfNormIn,bodyOut,bodyIn,sigma,yOrtho,yEff,intCells,intPoints" << endl;

    // Note (LK): this needs fixing 
    // loop over boundary cells
    //~ if (Pstream::nProcs() == 1)
    //~ {
        //~ List<List<intPoint>>& intInfos = lineIntInfoBoundary_->getIntPoints();

        //~ forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        //~ {
            //~ outFilePtr() << boundaryCells_[Pstream::myProcNo()][bCell].bCell_ << ","
                //~ << boundaryCells_[Pstream::myProcNo()][bCell].iCell_ << ","
                //~ << boundaryCells_[Pstream::myProcNo()][bCell].fCell1_ << ","
                //~ << boundaryCells_[Pstream::myProcNo()][bCell].fCell2_ << ","
                //~ << mesh_.C()[boundaryCells_[Pstream::myProcNo()][bCell].bCell_] << ","
                //~ << mesh_.C()[boundaryCells_[Pstream::myProcNo()][bCell].iCell_] << ","
                //~ << mesh_.C()[boundaryCells_[Pstream::myProcNo()][bCell].fCell1_] << ","
                //~ << mesh_.C()[boundaryCells_[Pstream::myProcNo()][bCell].fCell2_] << ","
                //~ << surfNorm_[boundaryCells_[Pstream::myProcNo()][bCell].bCell_] << ","
                //~ << surfNorm_[boundaryCells_[Pstream::myProcNo()][bCell].iCell_] << ","
                //~ << body_[boundaryCells_[Pstream::myProcNo()][bCell].bCell_] << ","
                //~ << body_[boundaryCells_[Pstream::myProcNo()][bCell].iCell_] << ","
                //~ << boundaryCells_[Pstream::myProcNo()][bCell].sigma_ << ","
                //~ << boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_ << ","
                //~ << boundaryCells_[Pstream::myProcNo()][bCell].yEff_ << ","
                //~ << intInfos[bCell][0].iCell_ << ","
                //~ << intInfos[bCell][0].iPoint_ << ","
                //~ << intInfos[bCell][1].iCell_ << ","
                //~ << intInfos[bCell][1].iPoint_ << ","
                //~ << intInfos[bCell][2].iCell_ << ","
                //~ << intInfos[bCell][2].iPoint_ << endl;
        //~ }
    //~ }

    // prepare file
    fileName = name + "_surface.dat";
    outFilePtr.reset(new OFstream(outDir/fileName));
    outFilePtr() << "cellI,cellCenter,surfNorm,body,sigma,order,intPoints,intCells" << endl;

    // loop over surface cells
    if (Pstream::nProcs() == 1)
    {
        List<List<intPoint>>& intInfos = lineIntInfoSurface_->getIntPoints();

        forAll(surfaceCells_[Pstream::myProcNo()], sCell)
        {
            outFilePtr () << surfaceCells_[Pstream::myProcNo()][sCell].sCell_ << ","
                << mesh_.C()[surfaceCells_[Pstream::myProcNo()][sCell].sCell_] << ","
                << surfNorm_[surfaceCells_[Pstream::myProcNo()][sCell].sCell_] << ","
                << body_[surfaceCells_[Pstream::myProcNo()][sCell].sCell_] << ","
                << surfaceCells_[Pstream::myProcNo()][sCell].sigma_ << ","
                << intInfos[sCell][0].iCell_ << ","
                << intInfos[sCell][0].iPoint_ << ","
                << intInfos[sCell][1].iCell_ << ","
                << intInfos[sCell][1].iPoint_ << ","
                << intInfos[sCell][2].iCell_ << ","
                << intInfos[sCell][2].iPoint_ << endl;
        }
    }
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
        label fProc = boundaryCells_[Pstream::myProcNo()][bCell].fProc1_;
        if (Pstream::myProcNo() != fProc)
        {
            fCellsToSync[fProc].append(boundaryCells_[Pstream::myProcNo()][bCell].fCell1_);
        }

        else
        {
            saveFreeCells.append(boundaryCells_[Pstream::myProcNo()][bCell].fCell1_);
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

// ************************************************************************* //
