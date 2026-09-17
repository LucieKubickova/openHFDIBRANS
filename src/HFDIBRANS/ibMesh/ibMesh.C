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
    openHFDIBRANS is licensed under the GNU LESSER GENERAL PUBLIC LICENSE
    (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the
    terms and conditions of version 3 of the GNU General Public License,
    supplemented by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIBRANS. If not, see
    <http://www.gnu.org/licenses/lgpl.html>.

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),
	Vít Večerník (2026-*)

\*---------------------------------------------------------------------------*/

#include "ibMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ibMesh::ibMesh
(
    const fvMesh& mesh,
    volScalarField& body
)
:
	mesh_(mesh),
	body_(body),
	yCorrected_(false),
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
	// Read HFDIBDEM dictionary
    stlName_ = HFDIBDEMDict_.lookupOrDefault<word>("stlName", "");
    cellSizeType_ =
		HFDIBDEMDict_.lookupOrDefault<word>("cellSizeType", "volumeRoot");
    valueL_ = HFDIBDEMDict_.lookupOrDefault<scalar>("sizeValue", 0.0);
    cutCellType_ =
		HFDIBDEMDict_.lookupOrDefault<word>("cutCellType", "cutCell");
	thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
	intSpan_ = readScalar(HFDIBDEMDict_.lookup("interfaceSpan"));
	sdBasedLambda_ =
		HFDIBDEMDict_.lookupOrDefault<bool>("sdBasedLamda", false);

	word geomModel =
		HFDIBDEMDict_.lookupOrDefault<word>("geomModel", "convex");
	bool genLambda =
		HFDIBDEMDict_.lookupOrDefault<bool>("generateLambda", false);

	if (!stlName_.empty())
	{
		// Read body stl
		bodySurfMesh_.reset
		(
			new triSurfaceMesh
			(
				IOobject
				(
					stlName_ + ".stl",
					mesh_.time().constant(),
					"triSurface",
					mesh_,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				)
			)
		);

		// Tri surface search
		triSurf_.reset(new triSurface(bodySurfMesh_()));
		triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));

		// Call to generate lambda
		initializeLambda(genLambda, geomModel);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ibMesh::~ibMesh()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void ibMesh::initializeLambda
(
	bool genLambda,
	word geomModel
)
{
	if (max(body_).value() < SMALL && genLambda)
	{
		Info<< "No initial lambda field found. Generating based on body: "
			<< stlName_ << endl;
		if (geomModel == "convex")
		{
			convexBody model
			(
				mesh_,
				thrSurf_,
				intSpan_,
				sdBasedLambda_,
				bodySurfMesh_,
				triSurfSearch_
			);
			model.generateLambda(body_);

			Info<< "Lambda field successfully generated" << endl;
		}
		else if (geomModel == "nonConvex")
		{
			nonConvexBody model
			(
				mesh_,
				thrSurf_,
				intSpan_,
				sdBasedLambda_,
				bodySurfMesh_,
				triSurfSearch_
			);
			model.generateLambda(body_);

			Info<< "Lambda field successfully generated" << endl;
		}
		else
		{
			FatalError
				<< "geomModel " << geomModel
				<<" not implemented" << exit(FatalError);
		}
	}
	else
	{
		Info<< "Initial lambda field provided" << endl;
	}

	// Update lambda values at the boundary
	Info<< "Correcting lambda boundary" << nl << endl;
	forAll(mesh_.boundaryMesh(), patchI)
	{
		const polyPatch& patch = mesh_.boundaryMesh()[patchI];

		if
		(
			!isA<emptyPolyPatch>(patch)
		 && !isA<processorPolyPatch>(patch)
		)
		{
			forAll(patch, faceI)
			{
				label owner = mesh_.faceOwner()[patch.start() + faceI];
				body_.boundaryFieldRef()[patchI][faceI] = body_[owner];
			}
		}
	}
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool ibMesh::pointInCell
(
    point pToCheck,
    label cToCheck
)
{
    const labelList& cellFaces = mesh_.cells()[cToCheck];
    forAll(cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        vector outNorm = mesh_.faceAreas()[fI];
        outNorm = (mesh_.faceOwner()[fI] == cToCheck) ? outNorm : (-1*outNorm);
        if (((pToCheck - mesh_.faceCentres()[fI]) & outNorm) > 0)
        {
            return false;
        }
    }
    return true;
}

// ------------------------------------------------------------------------- //

bool ibMesh::isWallCell
(
    label& cellI
)
{
    bool isWallCell = false;

    // Get wall patches
    DynamicList<label> wPatchIs;
    forAll(mesh_.boundary(), pI)
    {
        if (mesh_.boundary()[pI].type() == "wall")
        {
            wPatchIs.append(pI);
        }
    }

    // Loop over cell faces
    forAll(mesh_.cells()[cellI], f)
    {
        // Get face label
        label faceI = mesh_.cells()[cellI][f];

        if (faceI >= mesh_.owner().size())
        {
            bool wallFace = false;

            // Loop over patches of type wall
            forAll(wPatchIs, pI)
            {
                // Get patch label
                label patchI = wPatchIs[pI];

                // Get start and end face index
                label startI = mesh_.boundary()[patchI].start();
                label endI = startI + mesh_.boundary()[patchI].Cf().size();

                if (faceI >= startI && faceI < endI)
                {
                    wallFace = true;
                }
            }

            // Exclude wall faces
            if (wallFace)
            {
                isWallCell = true;
            }
        }
    }

    return isWallCell;
}

// ------------------------------------------------------------------------- //

bool ibMesh::isOnPatch
(
    label& cellI,
    word& patchName
)
{
    bool isOnPatch = false;

    // Get patch id
    const label patchI = mesh_.boundaryMesh().findPatchID(patchName);

    // Loop over cell faces
    forAll(mesh_.cells()[cellI], f)
    {
        // Get face label
        label faceI = mesh_.cells()[cellI][f];

        if (faceI >= mesh_.owner().size())
        {
            // Get start and end face index
            label startI = mesh_.boundary()[patchI].start();
            label endI = startI + mesh_.boundary()[patchI].Cf().size();

            if (faceI >= startI && faceI < endI)
            {
                isOnPatch = true;
            }
        }
    }

    return isOnPatch;
}

// ------------------------------------------------------------------------- //

label ibMesh::getFaceInDir
(
    label& cellI,
    vector& dir,
    label& prevFaceInDir
)
{
    // Prepare data
    label faceToReturn = -1;
    const labelList& cellFaces = mesh_.cells()[cellI];

    // Auxiliar scalar
    scalar dotProd = -GREAT;

    // Loop over cell faces
    forAll(cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        //~ vector outNorm = mesh_.Cf()[fI] - mesh_.C()[cellI];
        vector outNorm = mesh_.faceCentres()[fI] - mesh_.cellCentres()[cellI];
        outNorm /= mag(outNorm);

        //~ vector outNorm = (mesh_.faceOwner()[fI] == cellI)
            //~ ? mesh_.Sf()[fI] : (-1*mesh_.Sf()[fI]);
        //~ outNorm /= mag(outNorm); // LK: this should be there, no?

        scalar auxDotProd = outNorm & dir;
        if (auxDotProd > dotProd && fI != prevFaceInDir)
        {
            dotProd = auxDotProd;
            faceToReturn = fI;
        }
    }

    return faceToReturn;
}

// ------------------------------------------------------------------------- //

label ibMesh::getEdgeInDir
(
    label& faceI,
    label& cellI,
    vector& dir
)
{
    // Prepare data
    label edgeToReturn = -1;
    const labelList& faceEdges = mesh_.faceEdges()[faceI];

    // Auxiliar scalar
    scalar dotProd = -GREAT;

    // Loop over face edges
    forAll(faceEdges, edgeI)
    {
        // Get edge label
        label eI = faceEdges[edgeI];

        // Get edge nodes
        const label& own = mesh_.edges()[eI][0];
        const label& nei = mesh_.edges()[eI][1];

        vector Ce = 0.5*(mesh_.points()[own] + mesh_.points()[nei]);

        // Get direction to edge center
        vector outNorm = Ce - mesh_.C()[cellI];
        outNorm /= mag(outNorm);

        scalar auxDotProd = outNorm & dir;
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            edgeToReturn = eI;
        }
    }

    return edgeToReturn;
}

// ------------------------------------------------------------------------- //

label ibMesh::getVertInDir
(
    label& edgeI,
    label& cellI,
    vector& dir
)
{
    // Prepare data
    label vertexToReturn = -1;
    const edge& edgeVertices = mesh_.edges()[edgeI];

    // Auxiliar scalar
    scalar dotProd(-GREAT);

    // Loop over edge vertices
    forAll(edgeVertices, verI)
    {
        label vI = edgeVertices[verI];
        vector outNorm = mesh_.points()[vI] - mesh_.C()[cellI];
        outNorm /= mag(outNorm);

        scalar auxDotProd = outNorm & dir;
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            vertexToReturn = vI;
        }
    }

    return vertexToReturn;
}

// ------------------------------------------------------------------------- //

vector ibMesh::getClosestPoint
(
    vector ibPoint,
    intPoint& cPoint
)
{
    vector dir = cPoint.iPoint_ - ibPoint;
    dir /= mag(dir);

    vector dirToC = mesh_.C()[cPoint.iCell_] - ibPoint;

    return ibPoint + dir*(dirToC & dir);
}

// ------------------------------------------------------------------------- //

scalar ibMesh::createCutCellAndSurface
(
    label cellI,
    vector& normal,
    point& surfPoint
)
{
    scalar sArea = 0.0;

    if (cutCellType_ == "cutCell")
    {
        // Note (LK): original cut cell
        const cell& bCellSurf = mesh_.cells()[cellI];
        ibCutCell cCellSurf(mesh_, normal, surfPoint, bCellSurf);
        scalar yOrtho = cCellSurf.yOrtho();
		// Note (LK): creates the cut cell itself, should be as constructor

        // If the cell is uncut skip
        if (cCellSurf.faces().size() == 0)
        {
            Info<< "Warning: Uncut surface cell" << endl;
            return 0.0;
        }

        // Get area of cut face
        sArea = mag(cCellSurf.Sf()[cCellSurf.Sf().size() - 1]);
		// Note (LK): should be always the last one
    }

    // Note (LK): new cut cell, cutting edges by stl
    else if (cutCellType_ == "cutEdges")
    {
        // Prepare list of checked edges
        DynamicList<label> checkedEdges;

        // Get cell faces
        const labelList& cellFaces = mesh_.cells()[cellI];

        // Save points
        DynamicList<point> startPs;
        DynamicList<point> endPs;

        // Loop over cell faces
        forAll(cellFaces, fI)
        {
            // Get face label
            label faceI = cellFaces[fI];

            // Get face edges
            const labelList& faceEdges = mesh_.faceEdges()[faceI];

            // Loop over face edges
            forAll(faceEdges, eI)
            {
                // Look if already checked
                bool toInclude = true;
                forAll(checkedEdges, ceI)
                {
                    if (checkedEdges[ceI] == faceEdges[eI])
                    {
                        toInclude = false;
                        break;
                    }
                }

                // Break if already done
                if (!toInclude)
                {
                    continue;
                }
                else
                {
                    checkedEdges.append(faceEdges[eI]);
                }

                // Get edge
                const edge& e = mesh_.edges()[faceEdges[eI]];

                // Get points
                point sP = mesh_.points()[e.start()];
                point eP = mesh_.points()[e.end()];

                // Append
                startPs.append(sP);
                endPs.append(eP);
            }
        }

        // Prepare point fields
        pointField startPoints(startPs);
        pointField endPoints(endPs);

        // Try to find hit point with stl
        List<pointIndexHit> hitInfo;
        triSurfSearch_().findLine(startPoints, endPoints, hitInfo);

        // Get hit points
        DynamicList<point> cutPoints;
        forAll(hitInfo, hI)
        {
            if (hitInfo[hI].hit())
            {
                point cutPoint = hitInfo[hI].hitPoint();
                cutPoints.append(cutPoint);
            }
        }

        // Filter duplicate points
        DynamicList<point> uniquePoints;
        forAll(cutPoints, pI)
        {
            bool toAdd = true;

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

        // Calculate area
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
            Info<< "Warning: cell cut with " << uniquePoints.size()
				<< " points near " << mesh_.C()[cellI] << endl;
        }
    }

    else
    {
        FatalError
			<< "Surface area calculation type " << cutCellType_
			<< " not implemented" << exit(FatalError);
    }

    return sArea;
}

// ------------------------------------------------------------------------- //

scalar ibMesh::calculateTriangleArea
(
    point p0,
    point p1,
    point p2
)
{
    return mag(0.5*((p1 - p0)^(p2 - p1)));
}

// ------------------------------------------------------------------------- //

void ibMesh::createCutCellAndCenter
(
    label cellI,
    vector& surfNorm,
    point& surfPoint
)
{
    // Prepare list of checked edges
    DynamicList<label> checkedEdges;

    // Get cell faces
    const labelList& cellFaces = mesh_.cells()[cellI];

    // Save points
    DynamicList<point> startPs;
    DynamicList<point> endPs;

    // Loop over cell faces
    forAll(cellFaces, fI)
    {
        // Get face label
        label faceI = cellFaces[fI];

        // Get face edges
        const labelList& faceEdges = mesh_.faceEdges()[faceI];

        // Loop over face edges
        forAll(faceEdges, eI)
        {
            // Look if already checked
            bool toInclude = true;
            forAll(checkedEdges, ceI)
            {
                if (checkedEdges[ceI] == faceEdges[eI])
                {
                    toInclude = false;
                    break;
                }
            }

            // Break if already done
            if (!toInclude)
            {
                continue;
            }
            else
            {
                checkedEdges.append(faceEdges[eI]);
            }

            // Get edge
            const edge& e = mesh_.edges()[faceEdges[eI]];

            // Get points
            point sP = mesh_.points()[e.start()];
            point eP = mesh_.points()[e.end()];

            // Append
            startPs.append(sP);
            endPs.append(eP);
        }
    }

    // Prepare point fields
    pointField startPoints(startPs);
    pointField endPoints(endPs);

    // Try to find hit point with stl
    List<pointIndexHit> hitInfo;
    triSurfSearch_().findLine(startPoints, endPoints, hitInfo);

    // Get hit points
    DynamicList<point> cutPoints;
    forAll(hitInfo, hI)
    {
        if (hitInfo[hI].hit())
        {
            point cutPoint = hitInfo[hI].hitPoint();
            cutPoints.append(cutPoint);
        }
    }

    // Filter duplicate points
    DynamicList<point> uniquePoints;
    DynamicList<pointIndexHit> uniqueHitPoints;
    forAll(cutPoints, pI)
    {
        bool toAdd = true;

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

    // Get normals
    List<pointIndexHit> uniqueHitPointList(uniqueHitPoints);
    vectorField normalVectorField;

    // Get contact normal direction
    const triSurfaceMesh& ibTempMesh = bodySurfMesh_();
    ibTempMesh.getNormal(uniqueHitPointList, normalVectorField);

    // Calculate center and normal
    surfPoint *= 0.0;
    surfNorm *= 0.0;
    forAll(uniquePoints, uI)
    {
        surfPoint += uniquePoints[uI];
    }
    surfPoint /= uniquePoints.size();

    point closestPoint = vector::zero;
    scalar intDist = getCellSize(cellI, surfNorm);
    getClosestPointAndNormal(
        surfPoint,
        intDist*2*vector::one,
        closestPoint,
        surfNorm
    );

    return;
}

// ------------------------------------------------------------------------- //

void ibMesh::getClosestPointAndNormal
(
    const point& startPoint,
    const vector& span,
    point& closestPoint,
    vector& normal
)
{
    // Get nearest point on surface from contact center
    pointIndexHit ibPointIndexHit = triSurfSearch_().nearest(startPoint, span);
    List<pointIndexHit> ibPointIndexHitList(1, ibPointIndexHit);
    vectorField normalVectorField;

    // Get contact normal direction
    const triSurfaceMesh& ibTempMesh = bodySurfMesh_();
    ibTempMesh.getNormal(ibPointIndexHitList, normalVectorField);

    if(ibPointIndexHit.hit())
    {
        //~ normal = normalVectorField[0];
        closestPoint = ibPointIndexHit.hitPoint();
        normal = startPoint - closestPoint;
        normal /= mag(normal);
    }
    else
    {
        //~ FatalError
		//~		<< "Missing the closest point from " << startPoint
		//~		<< " to " << stlName_ << exit(FatalError);
        Info<< "Missing the closest point from " << startPoint
			<< " to " << stlName_ << endl;
    }
}

// ------------------------------------------------------------------------- //

scalar ibMesh::getCellSize
(
    label cellI,
    vector surfNorm
)
{
    scalar cellSize = 0.0;

    if (cellSizeType_ == "readSize")
    {
        cellSize = valueL_;
    }

    else if (cellSizeType_ == "volumeRoot" || mag(surfNorm) < SMALL)
    {
        cellSize = Foam::pow(mesh_.V()[cellI], 0.333);
    }

    else if (cellSizeType_ == "vertexBoundBox")
    {
        // Get cell vertices
        const labelList& cellVerts = mesh_.cellPoints()[cellI];

        // Prepare bounding box
        vector boundMin = mesh_.points()[cellVerts[0]];
        vector boundMax = mesh_.points()[cellVerts[0]];

        // Loop over vertices
        forAll(cellVerts, vI)
        {
            // Get vertex label
            label vertI = cellVerts[vI];

            // Get vertex point
            point vertP(mesh_.points()[vertI]);

            // Update bounding box
            boundMin.x() = min(boundMin.x(), vertP.x());
            boundMin.y() = min(boundMin.y(), vertP.y());
            boundMin.z() = min(boundMin.z(), vertP.z());

            boundMax.x() = max(boundMax.x(), vertP.x());
            boundMax.y() = max(boundMax.y(), vertP.y());
            boundMax.z() = max(boundMax.z(), vertP.z());
        }

        // Get bounding box size
        vector boundSize = boundMax - boundMin;

        // Get surface normal but in absolute values
        vector absSurfNorm =
			vector(mag(surfNorm.x()), mag(surfNorm.y()), mag(surfNorm.z()));

        // Get cell size in direction of surface normal
        cellSize = mag(boundSize & absSurfNorm);
    }


    return cellSize;
}

// ------------------------------------------------------------------------- //

void ibMesh::correctY
(
    volScalarField& y,
    bool recreate
)
{
    // If only lambda is provided, corrected in ibInterpolation
    if (sdBasedLambda_)
    {
        return;
    }

    // Skip if already done
    if (!recreate && yCorrected_)
    {
        return;
    }

    else
    {
        yCorrected_ = true;
    }

    // Set search distance span
    vector sDSpan = 4.0*(mesh_.bounds().max() - mesh_.bounds().min());

    // Loop over all cells
    forAll(y, cellI)
    {
        // Skip cells inside body
        if (body_[cellI] > 0.5)
        {
            continue;
        }

        // Prepare point and normal
        point closestPoint = vector::zero;
        vector surfNorm = vector::zero;

        // Get closest point and normal
        getClosestPointAndNormal(
            mesh_.C()[cellI],
            sDSpan,
            closestPoint,
            surfNorm
        );

        // Set corrected y
        y[cellI] = mag(mesh_.C()[cellI] - closestPoint);
    }
}


// ************************************************************************* //
