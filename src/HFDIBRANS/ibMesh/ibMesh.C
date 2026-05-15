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
    stlName_ = HFDIBDEMDict_.lookupOrDefault<word>("stlName", "");
    //~ readSurfNorm_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSurfaceNormal", false);
    //~ intSpan_ = readScalar(HFDIBDEMDict_.lookup("interfaceSpan"));
    //~ thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
    averageV_ = HFDIBDEMDict_.lookupOrDefault<bool>("averageVolume", false);
    readL_ = HFDIBDEMDict_.lookupOrDefault<bool>("readSize", false);
    valueL_ = HFDIBDEMDict_.lookupOrDefault<scalar>("sizeValue", 0.0); // LK: experimental
    sdBasedLambda_ = HFDIBDEMDict_.lookupOrDefault<bool>("sdBasedLambda", true);
    cutCellType_ = HFDIBDEMDict_.lookupOrDefault<word>("cutCellType", "cutCell");
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
bool ibMesh::isWallCell
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
bool ibMesh::isOnPatch
(
    label& cellI,
    word& patchName
)
{
    bool isOnPatch(false);

    // get patch id
    const label patchI = mesh_.boundaryMesh().findPatchID(patchName);

    // loop over cell faces
    forAll(mesh_.cells()[cellI], f)
    {
        // get face label
        label faceI = mesh_.cells()[cellI][f];

        if (faceI >= mesh_.owner().size())
        {
            // get start and end face index
            label startI = mesh_.boundary()[patchI].start();
            label endI = startI + mesh_.boundary()[patchI].Cf().size();

            if (faceI >= startI and faceI < endI)
            {
                isOnPatch = true;
            }
        }
    }

    return isOnPatch;
}

//---------------------------------------------------------------------------//
label ibMesh::getFaceInDir
(               
    label& cellI,
    vector& dir
)
{
    // prepare data
    label faceToReturn = -1;
    const labelList& cellFaces(mesh_.cells()[cellI]);

    // auxiliar scalar
    scalar dotProd(-GREAT);

    // loop over cell faces
    forAll(cellFaces, faceI)
    {    
        label fI = cellFaces[faceI];
        vector outNorm = mesh_.Cf()[fI] - mesh_.C()[cellI];
        outNorm /= mag(outNorm);
                    
        //~ vector outNorm = (mesh_.faceOwner()[fI] == cellI)
            //~ ? mesh_.Sf()[fI] : (-1*mesh_.Sf()[fI]);
        //~ outNorm /= mag(outNorm); // LK: this should be there, no?

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
label ibMesh::getEdgeInDir
(
    label& faceI,
    label& cellI,
    vector& dir
)
{
    // prepare data
    label edgeToReturn = -1;
    const labelList& faceEdges(mesh_.faceEdges()[faceI]);

    // auxiliar scalar
    scalar dotProd(-GREAT);

    // loop over face edges
    forAll(faceEdges, edgeI)
    {
        // get edge label
        label eI = faceEdges[edgeI];

        // get edge nodes
        const label& own = mesh_.edges()[eI][0];
        const label& nei = mesh_.edges()[eI][1];

        vector Ce = 0.5*(mesh_.points()[own] + mesh_.points()[nei]);

        // get direction to edge center
        vector outNorm = Ce - mesh_.C()[cellI];
        outNorm /= mag(outNorm);

        scalar auxDotProd(outNorm & dir);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            edgeToReturn = eI;
        }
    }

    return edgeToReturn;
}

//---------------------------------------------------------------------------//
label ibMesh::getVertInDir
(
    label& edgeI,
    label& cellI,
    vector& dir
)
{
    // prepare data
    label vertexToReturn = -1;
    const edge& edgeVertices(mesh_.edges()[edgeI]);

    // auxiliar scalar
    scalar dotProd(-GREAT);

    // loop over edge vertices
    forAll(edgeVertices, verI)
    {
        label vI = edgeVertices[verI];
        vector outNorm = mesh_.points()[vI] - mesh_.C()[cellI];
        outNorm /= mag(outNorm);

        scalar auxDotProd(outNorm & dir);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            vertexToReturn = vI;
        }
    }

    return vertexToReturn;
}

//---------------------------------------------------------------------------//
vector ibMesh::getClosestPoint
(
    vector ibPoint,
    intPoint& cPoint
)
{
    vector dir = cPoint.iPoint_ - ibPoint;
    dir /= mag(dir);

    vector dirToC = mesh_.C()[cPoint.iCell_] - ibPoint;

    return ibPoint + dir*(dirToC&dir);
}

//---------------------------------------------------------------------------//
scalar ibMesh::createCutCellAndSurface
(
    label cellI,
    vector& normal,
    point& surfPoint
)
{
    scalar sArea(0.0);

    if (cutCellType_ == "cutCell")
    {
        // Note (LK): original cut cell
        const cell& bCellSurf(mesh_.cells()[cellI]);
        ibCutCell cCellSurf(mesh_, normal, surfPoint, bCellSurf);
        scalar yOrtho = cCellSurf.yOrtho(); // Note (LK): creates the cut cell itself, should be as constructor

        // if the cell is uncut skip
        if (cCellSurf.faces().size() == 0)
        {
            Info << "Warning: Uncut surface cell" << endl;
            return 0.0;
        }

        // get area of cut face
        sArea = mag(cCellSurf.Sf()[cCellSurf.Sf().size()-1]); // Note (LK): should be always the last one
    }

    // Note (LK): new cut cell, cutting edges by stl
    else if (cutCellType_ == "cutEdges")
    {
        // prepare list of checked edges
        DynamicList<label> checkedEdges;

        // get cell faces
        const labelList& cellFaces(mesh_.cells()[cellI]);

        // save points
        DynamicList<point> startPs;
        DynamicList<point> endPs;

        // loop over cell faces
        forAll(cellFaces, fI)
        {
            // get face label
            label faceI = cellFaces[fI];

            // get face edges
            const labelList& faceEdges = mesh_.faceEdges()[faceI];

            // loop over face edges
            forAll(faceEdges, eI)
            {
                // look if already checked
                bool toInclude(true);
                forAll(checkedEdges, ceI)
                {
                    if (checkedEdges[ceI] == faceEdges[eI])
                    {
                        toInclude = false;
                        break;
                    }
                }

                // break if already don
                if (not toInclude)
                {
                    continue;
                }
                else
                {
                    checkedEdges.append(faceEdges[eI]);
                }

                // get edge
                const edge& e = mesh_.edges()[faceEdges[eI]];

                // get points
                point sP(mesh_.points()[e.start()]);
                point eP(mesh_.points()[e.end()]);

                // append
                startPs.append(sP);
                endPs.append(eP);
            }
        }

        // prepare point fields
        pointField startPoints(startPs);
        pointField endPoints(endPs);

        // try to find hit point with stl
        List<pointIndexHit> hitInfo;
        triSurfSearch_().findLine(startPoints, endPoints, hitInfo);

        // get hit points
        DynamicList<point> cutPoints;
        forAll(hitInfo, hI)
        {
            if (hitInfo[hI].hit())
            {
                point cutPoint = hitInfo[hI].hitPoint();
                cutPoints.append(cutPoint);
            }
        }

        // filter duplicate points
        DynamicList<point> uniquePoints;
        forAll(cutPoints, pI)
        {
            bool toAdd(true);

            forAll(uniquePoints, uI)
            {
                scalar dist = mag(cutPoints[pI] - uniquePoints[uI]);
                if (dist < SMALL)
                {
                    toAdd = false;
                }
            }

            if (toAdd)
            {
                uniquePoints.append(cutPoints[pI]);
            }
        }

        // calculate area
        if (uniquePoints.size() < 3)
        {
            sArea *= 0.0;
        }

        else if (uniquePoints.size() == 3)
        {
            sArea *= 0.0;

            point p0 = uniquePoints[0];
            point p1 = uniquePoints[1];
            point p2 = uniquePoints[2];

            sArea = calculateTriangleArea(p0, p1, p2);
        }

        else if (uniquePoints.size() == 4)
        {
            sArea *= 0.0;

            forAll(uniquePoints, pI)
            {
                point p0 = uniquePoints[pI % 4];
                point p1 = uniquePoints[(pI+1) % 4];
                point p2 = uniquePoints[(pI+2) % 4];
                sArea += calculateTriangleArea(p0, p1, p2);
            }

            sArea *= 0.5;
        }

        else
        {
            Info << "Warning: cell cut with " << uniquePoints.size() << " points near " << mesh_.C()[cellI] << endl;
        }
    }

    else
    {
        FatalError << "Surface area calculation type " << cutCellType_ << " not implemented" << exit(FatalError);
    }

    return sArea;
}

//---------------------------------------------------------------------------//
scalar ibMesh::calculateTriangleArea
(
    point p0,
    point p1,
    point p2
)
{
    return mag(0.5*((p1 - p0)^(p2 - p1)));
}

//---------------------------------------------------------------------------//
void ibMesh::createCutCellAndCenter
(
    label cellI,
    vector& surfNorm,
    point& surfPoint
)
{
    // prepare list of checked edges
    DynamicList<label> checkedEdges;

    // get cell faces
    const labelList& cellFaces(mesh_.cells()[cellI]);

    // save points
    DynamicList<point> startPs;
    DynamicList<point> endPs;

    // loop over cell faces
    forAll(cellFaces, fI)
    {
        // get face label
        label faceI = cellFaces[fI];

        // get face edges
        const labelList& faceEdges = mesh_.faceEdges()[faceI];

        // loop over face edges
        forAll(faceEdges, eI)
        {
            // look if already checked
            bool toInclude(true);
            forAll(checkedEdges, ceI)
            {
                if (checkedEdges[ceI] == faceEdges[eI])
                {
                    toInclude = false;
                    break;
                }
            }

            // break if already done
            if (not toInclude)
            {
                continue;
            }
            else
            {
                checkedEdges.append(faceEdges[eI]);
            }

            // get edge
            const edge& e = mesh_.edges()[faceEdges[eI]];

            // get points
            point sP(mesh_.points()[e.start()]);
            point eP(mesh_.points()[e.end()]);

            // append
            startPs.append(sP);
            endPs.append(eP);
        }
    }

    // prepare point fields
    pointField startPoints(startPs);
    pointField endPoints(endPs);

    // try to find hit point with stl
    List<pointIndexHit> hitInfo;
    triSurfSearch_().findLine(startPoints, endPoints, hitInfo);

    // get hit points
    DynamicList<point> cutPoints;
    forAll(hitInfo, hI)
    {
        if (hitInfo[hI].hit())
        {
            point cutPoint = hitInfo[hI].hitPoint();
            cutPoints.append(cutPoint);
        }
    }

    // filter duplicate points
    DynamicList<point> uniquePoints;
    DynamicList<pointIndexHit> uniqueHitPoints;
    forAll(cutPoints, pI)
    {
        bool toAdd(true);

        forAll(uniquePoints, uI)
        {
            scalar dist = mag(cutPoints[pI] - uniquePoints[uI]);
            if (dist < SMALL)
            {
                toAdd = false;
            }
        }

        if (toAdd)
        {
            uniquePoints.append(cutPoints[pI]);
            uniqueHitPoints.append(hitInfo[pI]);
        }
    }

    // get normals
    List<pointIndexHit> uniqueHitPointList(uniqueHitPoints);
    vectorField normalVectorField;

    // get contact normal direction
    const triSurfaceMesh& ibTempMesh(bodySurfMesh_());
    ibTempMesh.getNormal(uniqueHitPointList,normalVectorField);

    // calculate center and normal
    surfPoint *= 0.0;
    surfNorm *= 0.0;
    forAll(uniquePoints, uI)
    {
        surfPoint += uniquePoints[uI];
        //~ surfNorm += normalVectorField[uI];
    }
    surfPoint /= uniquePoints.size();

    point closestPoint(vector::zero);
    scalar intDist = Foam::pow(mesh_.V()[cellI],0.333);
    getClosestPointAndNormal(
        surfPoint,
        intDist*2*vector::one,
        closestPoint,
        surfNorm
    );

    //~ surfNorm /= uniquePoints.size();

    return;
}

//---------------------------------------------------------------------------//
void ibMesh::getClosestPointAndNormal
(
    const point& startPoint,
    const vector& span,
    point& closestPoint,
    vector& normal
)
{
    // get nearest point on surface from contact center
    pointIndexHit ibPointIndexHit = triSurfSearch_().nearest(startPoint, span);
    List<pointIndexHit> ibPointIndexHitList(1,ibPointIndexHit);
    vectorField normalVectorField;

    // get contact normal direction
    const triSurfaceMesh& ibTempMesh(bodySurfMesh_());
    ibTempMesh.getNormal(ibPointIndexHitList,normalVectorField);

    if(ibPointIndexHit.hit())
    {
        //~ normal = normalVectorField[0];
        closestPoint = ibPointIndexHit.hitPoint();
        normal = startPoint - closestPoint;
        normal /= mag(normal);
    }
    else
    {
        FatalError << "Missing the closest point from " << startPoint << " to " << stlName_ << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
scalar ibMesh::getCellSize
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

// ************************************************************************* //
