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

#include "ibMesh.H"

using namespace Foam;

//---------------------------------------------------------------------------//
ibMesh::ibMesh
(
    const fvMesh& mesh,
    const volScalarField& body
)
:
mesh_(mesh),
body_(body),
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
    //~ boundarySearch_ = HFDIBDEMDict_.lookupOrDefault<word>("boundarySearch", "face");
    stlName_ = HFDIBDEMDict_.lookupOrDefault<word>("stlName", "");
    //~ excludeWalls_ = HFDIBDEMDict_.lookupOrDefault<bool>("excludeWalls", false);
    //~ excludePatch_ = HFDIBDEMDict_.lookupOrDefault<word>("excludePatch", "none");
    //~ readSurfNorm_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSurfaceNormal", false);
    //~ intSpan_ = readScalar(HFDIBDEMDict_.lookup("interfaceSpan"));
    //~ thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
    //~ averageV_ = HFDIBDEMDict_.lookupOrDefault<bool>("averageVolume", false);
    //~ readL_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSize", false);
    //~ valueL_ = HFDIBDEMDict_.lookupOrDefault<scalar>("sizeValue", 0.0); // LK: experimental
    sdBasedLambda_ = HFDIBDEMDict_.lookupOrDefault<bool>("sdBasedLambda", true);
    //~ surfAreaType_ = HFDIBDEMDict_.lookupOrDefault<word>("surfAreaType", "cutCell");
    //~ yFromCutEdges_ = HFDIBDEMDict_.lookupOrDefault<bool>("yFromCutEdges", false); // LK: experimental

    if (!sdBasedLambda_)
    {
        stlPath_ = "constant/triSurface/" + stlName_ + ".stl";

        // read stl
        bodySurfMesh_.reset(new triSurfaceMesh
        (
            IOobject
            (
                stlPath_,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ));

        // tri surface search
        triSurf_.reset(new triSurface(bodySurfMesh_()));
        triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
    }

    //~ // compute average cell volume
    //~ VAve_ = 0.0;
    //~ if (averageV_)
    //~ {
        //~ forAll(mesh_.V(), i)
        //~ {
            //~ VAve_ += mesh_.V()[i];
        //~ }
        //~ VAve_ /= mesh_.V().size();
    //~ }

    // calculate surface normals
    //~ calculateSurfNorm();
}

//---------------------------------------------------------------------------//
ibMesh::~ibMesh()
{
}

//---------------------------------------------------------------------------//
bool ibMesh::pointInCell
(           
    point pToCheck,
    label cToCheck
)
{           
    const labelList& cellFaces(mesh_.cells()[cToCheck]);
    forAll(cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        vector outNorm = mesh_.Sf()[fI];
        outNorm = (mesh_.faceOwner()[fI] == cToCheck) ? outNorm : (-1*outNorm);
                
        if (((pToCheck - mesh_.Cf()[fI]) & outNorm) > 0)
        {
            return false;
        }
    }       
    return true;
}

//---------------------------------------------------------------------------//
//~ void ibMesh::findNeighborInBody
//~ (
    //~ label& cellI,
    //~ scalar threshold,
    //~ label& iCell,
    //~ label& iFace,
    //~ label& iProc,
    //~ bool& isBoundary
//~ )
//~ {
    //~ // search
    //~ if (boundarySearch_ == "vertex")
    //~ {
        //~ if (Pstream::nProcs() == 1) // Note (LK): not really working single core version, gonna be removed
        //~ {
            //~ forAll(mesh_.cellPoints()[cellI], pID) // vertex neighbours
            //~ {
                //~ label pointI = mesh_.cellPoints()[cellI][pID];
    
                //~ forAll(mesh_.pointCells()[pointI], cI)
                //~ {
                    //~ if (body_[mesh_.pointCells()[pointI][cI]] >= threshold)
                    //~ {
                        //~ isBoundary = true;
                        //~ iCell = mesh_.pointCells()[pointI][cI];
                        //~ break;
                    //~ }
                //~ }
            //~ }
        //~ }

        //~ else
        //~ {
            //~ // Note (LK): should be fixed later
            //~ FatalError << "Boundary cell search " << boundarySearch_ << " not implemented in parallel" << exit(FatalError);
        //~ }
    //~ }

    //~ else if (boundarySearch_ == "edge")
    //~ {
        //~ // get the best face, edge and vertex
        //~ label faceI = getFaceInDir(cellI, -1*surfNorm_[cellI]);
        //~ label edgeI = getEdgeInDir(faceI, cellI, -1*surfNorm_[cellI]);
        //~ //~ label vertI = getVertInDir(edgeI, cellI, -1*surfNorm_[cellI]);

        //~ label faceNI(-1);
        //~ label edgeNI(-1);
    
        //~ // check for non-internal cells
        //~ if (!mesh_.isInternalFace(faceI))
        //~ {
            //~ //~ // get the patch the face belongs to
            //~ //~ label facePatchI(mesh_.boundaryMesh().whichPatch(faceI));
            //~ //~ const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchI];
    
            //~ //~ // check if it is a processor boundary
            //~ //~ if (cPatch.type() == "processor")
            //~ //~ {
                //~ //~ // get the processor patch
                //~ //~ const processorPolyPatch& procPatch
                    //~ //~ = refCast<const processorPolyPatch>(cPatch);
    
                //~ //~ // get the neighboring processor id
                //~ //~ label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                    //~ //~ ? procPatch.neighbProcNo() : procPatch.myProcNo();
    
                //~ //~ // access neighbor processor patch
                //~ //~ const scalarField& bodyN = body_.boundaryField()[facePatchI].patchNeighbourField();
    
                //~ //~ // get local face value
                //~ //~ label localI = cPatch.whichFace(faceI);
    
                //~ //~ // acces neighbor value
                //~ //~ if (bodyN[localI] >= threshold)
                //~ //~ {
                    //~ //~ isBoundary = true;
                    //~ //~ iFace = localI;
                    //~ //~ iProc = sProc;
                //~ //~ }
            //~ //~ }
        //~ }
    
        //~ else
        //~ {
            //~ // get the face neighbor
            //~ label owner(mesh_.owner()[faceI]);
            //~ label neighbor(mesh_.neighbour()[faceI]);
            //~ faceNI = (cellI == owner) ? neighbor : owner;

            //~ // get the edge neighbor
            //~ bool found(false);
            //~ forAll(mesh_.edgeCells()[edgeI], cI)
            //~ {
                //~ // get cell label
                //~ edgeNI = mesh_.edgeCells()[edgeI][cI];

                //~ // skip current cell and the face neighbor
                //~ if (edgeNI == cellI or edgeNI == faceNI)
                //~ {
                    //~ continue;
                //~ }

                //~ // check if vertex neighbor is neighbor with face neighbor O:)
                //~ forAll(mesh_.cellCells()[faceNI], pI)
                //~ {
                    //~ label posNI = mesh_.cellCells()[faceNI][pI];
                    //~ if (posNI == edgeNI)
                    //~ {
                        //~ found = true;
                    //~ }
                //~ }

                //~ if (found)
                //~ {
                    //~ break;
                //~ }
            //~ }
    
            //~ // check lambda field
            //~ // LK: needs check if edgeNI found
            //~ // LK: needs check which more in dir
            //~ if (body_[faceNI] >= threshold)
            //~ {
                //~ isBoundary = true;
                //~ iCell = faceNI;
                //~ iProc = Pstream::myProcNo();
            //~ }

            //~ else if (body_[edgeNI] >= threshold)
            //~ {
                //~ isBoundary = true;
                //~ iCell = edgeNI;
                //~ iProc = Pstream::myProcNo();
            //~ }
        //~ }

        //~ // Note (LK): should be fixed later
        //~ FatalError << "Boundary cell search " << boundarySearch_ << " not finished" << exit(FatalError);
    //~ }
    
    //~ else if (boundarySearch_ == "face")
    //~ {
        //~ // get labels
        //~ // Note (LK): surf norm should be from STL file if stated
        //~ label faceI = getFaceInDir(cellI, -1*surfNorm_[cellI]);
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
    
                //~ // access neighbor processor patch
                //~ const tmp<scalarField> tbodyN(body_.boundaryField()[facePatchI].patchNeighbourField());
                //~ const scalarField& bodyN = tbodyN();
    
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
        //~ }
    
        //~ else
        //~ {
            //~ // get the neighbor cell label
            //~ label owner(mesh_.owner()[faceI]);
            //~ label neighbor(mesh_.neighbour()[faceI]);
            //~ nI = (cellI == owner) ? neighbor : owner;

            //~ // check lambda field
            //~ if (body_[nI] >= threshold)
            //~ {
                //~ isBoundary = true;
                //~ iCell = nI;
                //~ iProc = Pstream::myProcNo();
            //~ }
        //~ }
    //~ }
    
    //~ else
    //~ {
        //~ FatalError << "Boundary cell search " << boundarySearch_ << " not implemented" << exit(FatalError);
    //~ }
//~ }

//~ //---------------------------------------------------------------------------//
//~ label ibMesh::getFaceInDir
//~ (
    //~ label cellI,
    //~ vector dir
//~ )
//~ {
    //~ // prepare data
    //~ label faceToReturn = -1;
    //~ const labelList& cellFaces(mesh_.cells()[cellI]);

    //~ // auxiliar scalar
    //~ scalar dotProd(-GREAT);

    //~ // loop over cell faces
    //~ forAll(cellFaces, faceI)
    //~ {
        //~ label fI = cellFaces[faceI];
        //~ vector outNorm = mesh_.Cf()[fI] - mesh_.C()[cellI];
        //~ outNorm /= mag(outNorm);

        //~ //~ vector outNorm = (mesh_.faceOwner()[fI] == cellI)
            //~ //~ ? mesh_.Sf()[fI] : (-1*mesh_.Sf()[fI]);
        //~ //~ outNorm /= mag(outNorm); // LK: this should be there, no?

        //~ scalar auxDotProd(outNorm & dir);
        //~ if (auxDotProd > dotProd)
        //~ {
            //~ dotProd = auxDotProd;
            //~ faceToReturn = fI;
        //~ }
    //~ }

    //~ return faceToReturn;
//~ }

//~ //---------------------------------------------------------------------------//
//~ label ibMesh::getEdgeInDir
//~ (
    //~ label faceI,
    //~ label cellI,
    //~ vector dir
//~ )
//~ {
    //~ // prepare data
    //~ label edgeToReturn = -1;
    //~ const labelList& faceEdges(mesh_.faceEdges()[faceI]);

    //~ // auxiliar scalar
    //~ scalar dotProd(-GREAT);

    //~ // loop over face edges
    //~ forAll(faceEdges, edgeI)
    //~ {
        //~ // get edge label
        //~ label eI = faceEdges[edgeI];

        //~ // get edge nodes
        //~ const label& own = mesh_.edges()[eI][0];
        //~ const label& nei = mesh_.edges()[eI][1];

        //~ vector Ce = 0.5*(mesh_.points()[own] + mesh_.points()[nei]);

        //~ // get direction to edge center
        //~ vector outNorm = Ce - mesh_.C()[cellI];
        //~ outNorm /= mag(outNorm);

        //~ scalar auxDotProd(outNorm & dir);
        //~ if (auxDotProd > dotProd)
        //~ {
            //~ dotProd = auxDotProd;
            //~ edgeToReturn = eI;
        //~ }
    //~ }

    //~ return edgeToReturn;
//~ }

//~ //---------------------------------------------------------------------------//
//~ label ibMesh::getVertInDir
//~ (
    //~ label edgeI,
    //~ label cellI,
    //~ vector dir
//~ )
//~ {
    //~ // prepare data
    //~ label vertexToReturn = -1;
    //~ const edge& edgeVertices(mesh_.edges()[edgeI]);

    //~ // auxiliar scalar
    //~ scalar dotProd(-GREAT);

    //~ // loop over edge vertices
    //~ forAll(edgeVertices, verI)
    //~ {
        //~ label vI = edgeVertices[verI];
        //~ vector outNorm = mesh_.points()[vI] - mesh_.C()[cellI];
        //~ outNorm /= mag(outNorm);

        //~ scalar auxDotProd(outNorm & dir);
        //~ if (auxDotProd > dotProd)
        //~ {
            //~ dotProd = auxDotProd;
            //~ vertexToReturn = vI;
        //~ }
    //~ }

    //~ return vertexToReturn;
//~ }

//~ //---------------------------------------------------------------------------//
//~ bool ibMesh::isWallCell
//~ (
    //~ label& cellI
//~ )
//~ {
    //~ bool isWallCell(true);

    //~ // get wall patches
    //~ DynamicList<label> wPatchIs;
    //~ forAll(mesh_.boundary(), pI)
    //~ {
        //~ if (mesh_.boundary()[pI].type() == "wall")
        //~ {
            //~ wPatchIs.append(pI);
        //~ }
    //~ }

    //~ // loop over cell faces
    //~ forAll(mesh_.cells()[cellI], f)
    //~ {
        //~ // get face label
        //~ label faceI = mesh_.cells()[cellI][f];
    
        //~ if (faceI >= mesh_.owner().size())
        //~ {
            //~ bool wallFace(false);
    
            //~ // loop over patches of type wall
            //~ forAll(wPatchIs, pI)
            //~ {
                //~ // get patch label
                //~ label patchI = wPatchIs[pI];
    
                //~ // get start and end face index
                //~ label startI = mesh_.boundary()[patchI].start();
                //~ label endI = startI + mesh_.boundary()[patchI].Cf().size();
    
                //~ if (faceI >= startI and faceI < endI)
                //~ {
                    //~ wallFace = true;
                //~ }
            //~ }
    
            //~ // exclude wall faces
            //~ if (wallFace)
            //~ {
                //~ isWallCell = false;
            //~ }
        //~ }
    //~ }

    //~ return isWallCell;
//~ }

//~ //---------------------------------------------------------------------------//
//~ bool ibMesh::isOnPatch
//~ (
    //~ label& cellI
//~ )
//~ {
    //~ bool isOnPatch(false);

    //~ // get patch id
    //~ const label patchI = mesh_.boundaryMesh().findPatchID(excludePatch_);

    //~ // loop over cell faces
    //~ forAll(mesh_.cells()[cellI], f)
    //~ {
        //~ // get face label
        //~ label faceI = mesh_.cells()[cellI][f];
    
        //~ if (faceI >= mesh_.owner().size())
        //~ {
            //~ // get start and end face index
            //~ label startI = mesh_.boundary()[patchI].start();
            //~ label endI = startI + mesh_.boundary()[patchI].Cf().size();
    
            //~ if (faceI >= startI and faceI < endI)
            //~ {
                //~ isOnPatch = true;
            //~ }
        //~ }
    //~ }

    //~ return isOnPatch;
//~ }

//~ //---------------------------------------------------------------------------//
//~ scalar ibMesh::createCutCellAndSurface
//~ (
    //~ label cellI,
    //~ vector& normal,
    //~ point& surfPoint
//~ )
//~ {
    //~ scalar sArea(0.0);

    //~ if (surfAreaType_ == "cutCell")
    //~ {
        //~ // Note (LK): original cut cell
        //~ const cell& bCellSurf(mesh_.cells()[cellI]);
        //~ ibCutCell cCellSurf(mesh_, normal, surfPoint, bCellSurf);
        //~ scalar yOrtho = cCellSurf.yOrtho(); // Note (LK): creates the cut cell itself, should be as constructor
        
        //~ // if the cell is uncut skip
        //~ if (cCellSurf.faces().size() == 0)
        //~ {
            //~ Info << "Warning: Uncut surface cell" << endl;
            //~ return 0.0;
        //~ }
        
        //~ // get area of cut face
        //~ sArea = mag(cCellSurf.Sf()[cCellSurf.Sf().size()-1]); // Note (LK): should be always the last one
    //~ }

    //~ // Note (LK): new cut cell, cutting edges by stl
    //~ else if (surfAreaType_ == "cutEdges")
    //~ {
        //~ // prepare list of checked edges
        //~ DynamicList<label> checkedEdges;
        
        //~ // get cell faces
        //~ const labelList& cellFaces(mesh_.cells()[cellI]);

        //~ // save points
        //~ DynamicList<point> startPs;
        //~ DynamicList<point> endPs;

        //~ // loop over cell faces
        //~ forAll(cellFaces, fI)
        //~ {
            //~ // get face label
            //~ label faceI = cellFaces[fI];

            //~ // get face edges
            //~ const labelList& faceEdges = mesh_.faceEdges()[faceI];

            //~ // loop over face edges
            //~ forAll(faceEdges, eI)
            //~ {
                //~ // look if already checked
                //~ bool toInclude(true);
                //~ forAll(checkedEdges, ceI)
                //~ {
                    //~ if (checkedEdges[ceI] == faceEdges[eI])
                    //~ {
                        //~ toInclude = false;
                        //~ break;
                    //~ }
                //~ }

                //~ // break if already don
                //~ if (not toInclude)
                //~ {
                    //~ continue;
                //~ }
                //~ else
                //~ {
                    //~ checkedEdges.append(faceEdges[eI]);
                //~ }

                //~ // get edge
                //~ const edge& e = mesh_.edges()[faceEdges[eI]];

                //~ // get points
                //~ point sP(mesh_.points()[e.start()]);
                //~ point eP(mesh_.points()[e.end()]);

                //~ // append
                //~ startPs.append(sP);
                //~ endPs.append(eP);
            //~ }
        //~ }

        //~ // prepare point fields
        //~ pointField startPoints(startPs);
        //~ pointField endPoints(endPs);

        //~ // try to find hit point with stl
        //~ List<pointIndexHit> hitInfo;
        //~ triSurfSearch_().findLine(startPoints, endPoints, hitInfo);

        //~ // get hit points
        //~ DynamicList<point> cutPoints;
        //~ forAll(hitInfo, hI)
        //~ {
            //~ if (hitInfo[hI].hit())
            //~ {
                //~ point cutPoint = hitInfo[hI].hitPoint();
                //~ cutPoints.append(cutPoint);
            //~ }
        //~ }

        //~ // filter duplicate points
        //~ DynamicList<point> uniquePoints;
        //~ forAll(cutPoints, pI)
        //~ {
            //~ bool toAdd(true);

            //~ forAll(uniquePoints, uI)
            //~ {
                //~ scalar dist = mag(cutPoints[pI] - uniquePoints[uI]);
                //~ if (dist < SMALL)
                //~ {
                    //~ toAdd = false;
                //~ }
            //~ }

            //~ if (toAdd)
            //~ {
                //~ uniquePoints.append(cutPoints[pI]);
            //~ }
        //~ }

        //~ // calculate area
        //~ if (uniquePoints.size() < 3)
        //~ {
            //~ sArea *= 0.0;
        //~ }

        //~ else if (uniquePoints.size() == 3)
        //~ {
            //~ sArea *= 0.0;

            //~ point p0 = uniquePoints[0];
            //~ point p1 = uniquePoints[1];
            //~ point p2 = uniquePoints[2];

            //~ sArea = calculateTriangleArea(p0, p1, p2);
        //~ }

        //~ else if (uniquePoints.size() == 4)
        //~ {
            //~ sArea *= 0.0;

            //~ forAll(uniquePoints, pI)
            //~ {
                //~ point p0 = uniquePoints[pI % 4];
                //~ point p1 = uniquePoints[(pI+1) % 4];
                //~ point p2 = uniquePoints[(pI+2) % 4];
                //~ sArea += calculateTriangleArea(p0, p1, p2);
            //~ }

            //~ sArea *= 0.5;
        //~ }

        //~ else
        //~ {
            //~ Info << "Warning: cell cut with " << uniquePoints.size() << " points near " << mesh_.C()[cellI] << endl;
        //~ }
    //~ }

    //~ else
    //~ {
        //~ FatalError << "Surface area calculation type " << surfAreaType_ << " not implemented" << exit(FatalError);
    //~ }

    //~ return sArea;
//~ }

//~ //---------------------------------------------------------------------------//
//~ scalar ibMesh::calculateTriangleArea
//~ (
    //~ point p0,
    //~ point p1,
    //~ point p2
//~ )
//~ {
    //~ return mag(0.5*((p1 - p0)^(p2 - p1)));
//~ }

//~ //---------------------------------------------------------------------------//
//~ void ibMesh::createCutCellAndCenter
//~ (
    //~ label cellI,
    //~ vector& surfNorm,
    //~ point& surfPoint
//~ )
//~ {
    //~ // prepare list of checked edges
    //~ DynamicList<label> checkedEdges;
    
    //~ // get cell faces
    //~ const labelList& cellFaces(mesh_.cells()[cellI]);

    //~ // save points
    //~ DynamicList<point> startPs;
    //~ DynamicList<point> endPs;

    //~ // loop over cell faces
    //~ forAll(cellFaces, fI)
    //~ {
        //~ // get face label
        //~ label faceI = cellFaces[fI];

        //~ // get face edges
        //~ const labelList& faceEdges = mesh_.faceEdges()[faceI];

        //~ // loop over face edges
        //~ forAll(faceEdges, eI)
        //~ {
            //~ // look if already checked
            //~ bool toInclude(true);
            //~ forAll(checkedEdges, ceI)
            //~ {
                //~ if (checkedEdges[ceI] == faceEdges[eI])
                //~ {
                    //~ toInclude = false;
                    //~ break;
                //~ }
            //~ }

            //~ // break if already done
            //~ if (not toInclude)
            //~ {
                //~ continue;
            //~ }
            //~ else
            //~ {
                //~ checkedEdges.append(faceEdges[eI]);
            //~ }

            //~ // get edge
            //~ const edge& e = mesh_.edges()[faceEdges[eI]];

            //~ // get points
            //~ point sP(mesh_.points()[e.start()]);
            //~ point eP(mesh_.points()[e.end()]);

            //~ // append
            //~ startPs.append(sP);
            //~ endPs.append(eP);
        //~ }
    //~ }

    //~ // prepare point fields
    //~ pointField startPoints(startPs);
    //~ pointField endPoints(endPs);

    //~ // try to find hit point with stl
    //~ List<pointIndexHit> hitInfo;
    //~ triSurfSearch_().findLine(startPoints, endPoints, hitInfo);

    //~ // get hit points
    //~ DynamicList<point> cutPoints;
    //~ forAll(hitInfo, hI)
    //~ {
        //~ if (hitInfo[hI].hit())
        //~ {
            //~ point cutPoint = hitInfo[hI].hitPoint();
            //~ cutPoints.append(cutPoint);
        //~ }
    //~ }

    //~ // filter duplicate points
    //~ DynamicList<point> uniquePoints;
    //~ DynamicList<pointIndexHit> uniqueHitPoints;
    //~ forAll(cutPoints, pI)
    //~ {
        //~ bool toAdd(true);

        //~ forAll(uniquePoints, uI)
        //~ {
            //~ scalar dist = mag(cutPoints[pI] - uniquePoints[uI]);
            //~ if (dist < SMALL)
            //~ {
                //~ toAdd = false;
            //~ }
        //~ }

        //~ if (toAdd)
        //~ {
            //~ uniquePoints.append(cutPoints[pI]);
            //~ uniqueHitPoints.append(hitInfo[pI]);
        //~ }
    //~ }

    //~ // get normals
    //~ List<pointIndexHit> uniqueHitPointList(uniqueHitPoints);
    //~ vectorField normalVectorField;

    //~ // get contact normal direction
    //~ const triSurfaceMesh& ibTempMesh(bodySurfMesh_());
    //~ ibTempMesh.getNormal(uniqueHitPointList,normalVectorField);

    //~ // calculate center and normal
    //~ surfPoint *= 0.0;
    //~ surfNorm *= 0.0;
    //~ forAll(uniquePoints, uI)
    //~ {
        //~ surfPoint += uniquePoints[uI];
        //~ //~ surfNorm += normalVectorField[uI];
    //~ }
    //~ surfPoint /= uniquePoints.size();

    //~ point closestPoint(vector::zero);
    //~ scalar intDist = Foam::pow(mesh_.V()[cellI],0.333);
    //~ getClosestPointAndNormal(
        //~ surfPoint,
        //~ intDist*2*vector::one,
        //~ closestPoint,
        //~ surfNorm 
    //~ );

    //~ //~ surfNorm /= uniquePoints.size();

    //~ return;
//~ }

//~ //---------------------------------------------------------------------------//
//~ scalar ibMesh::getCellSize
//~ (
    //~ label cellI
//~ )
//~ {
    //~ scalar cellSize(0.0);

    //~ if (averageV_)
    //~ {
        //~ cellSize = Foam::pow(VAve_, 0.333);
    //~ }
    //~ else if (readL_)
    //~ {
        //~ cellSize = valueL_;
    //~ }
    //~ else
    //~ {
        //~ cellSize = Foam::pow(mesh_.V()[cellI], 0.333);
    //~ }

    //~ return cellSize;
//~ }

//~ //---------------------------------------------------------------------------//
//~ void ibMesh::getClosestPointAndNormal
//~ (
    //~ const point& startPoint,
    //~ const vector& span,
    //~ point& closestPoint,
    //~ vector& normal
//~ )
//~ {
    //~ // get nearest point on surface from contact center
    //~ pointIndexHit ibPointIndexHit = triSurfSearch_().nearest(startPoint, span);
    //~ List<pointIndexHit> ibPointIndexHitList(1,ibPointIndexHit);
    //~ vectorField normalVectorField;

    //~ // get contact normal direction
    //~ const triSurfaceMesh& ibTempMesh(bodySurfMesh_());
    //~ ibTempMesh.getNormal(ibPointIndexHitList,normalVectorField);

    //~ if(ibPointIndexHit.hit())
    //~ {
        //~ //~ normal = normalVectorField[0];
        //~ closestPoint = ibPointIndexHit.hitPoint();
        //~ normal = startPoint - closestPoint;
        //~ normal /= mag(normal);
    //~ }
    //~ else
    //~ {
        //~ FatalError << "Missing the closest point from " << startPoint << " to " << stlName_ << exit(FatalError);
    //~ }
//~ }

// ************************************************************************* //
