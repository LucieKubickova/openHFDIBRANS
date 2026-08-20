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

#include "HashSet.H"
#include "HashTable.H"
#include "Pstream.H"
#include "error.H"
#include "pointIndexHit.H"
#include "vector.H"

#include "stlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stlModel::stlModel
(
	const fvMesh& mesh,
	scalar thrSurf,
	scalar intSpan,
	bool sdBasedLambda,
	autoPtr<triSurfaceMesh>& triSurfMesh,
	autoPtr<triSurfaceSearch>& triSurfSearch
)
:
	mesh_(mesh),
	thrSurf_(thrSurf),
	intSpan_(intSpan),
	sdBasedLambda_(sdBasedLambda),
	triSurfMesh_(triSurfMesh),
	triSurfSearch_(triSurfSearch),
	geometricD_(mesh_.geometricD()),
	meshBounds_(mesh_.points(), false),
	cellToStart_(0),
	cellPoints_(mesh_.nCells())
{
	internalCells_.setSize(Pstream::nProcs());

	forAll(mesh.cells(), cellI)
	{
		const labelList& ptIndices = mesh.cellPoints()[cellI];
		cellPoints_[cellI].setSize(ptIndices.size());
		forAll(ptIndices, pI)
		{
			cellPoints_[cellI][pI] = mesh.points()[ptIndices[pI]];
		}
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

stlModel::~stlModel()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

labelList stlModel::classifyCell
(
	label cellI,
	bool& centreInside,
	boolList& verticesInside,
	volScalarField& lambda
)
{
	// Vertex weight contribution:
	// each vertex = 0.5/N, centre = 0.5 -> 1.0 in total
	const scalar vertexWeight = 0.5/verticesInside.size();
	scalar cBody = 0;

	forAll(verticesInside, vI)
	{
		if (verticesInside[vI])
		{
			cBody += vertexWeight;
		}
	}
	if (centreInside)
	{
		cBody += 0.5;
	}

	labelList neighbours;

	if (cBody > thrSurf_)
	{
		if (cBody > (1.0 - thrSurf_))
		{
			// Fully internal cell
			internalCells_[Pstream::myProcNo()].append(cellI);
		}
		else if (sdBasedLambda_)
		{
			// Replace vertex weight estimate with tanh profiling
			// based on signed distance to nearest STL surface
			const vector sDSpan
			(
				4.0*(mesh_.bounds().max() - mesh_.bounds().min())
			);
			pointIndexHit hit
			(
				triSurfSearch_().nearest(mesh_.C()[cellI], sDSpan)
			);

			scalar dist = 0;
			if (hit.hit())
			{
				dist = mag(hit.hitPoint() - mesh_.C()[cellI]);
			}
			else
			{
				Info << "Missed the closest point from cell "
					 << cellI << " to STL body"  << endl;
			}

			const scalar cellSize = Foam::pow(mesh_.V()[cellI], 0.333);

			if (centreInside)
			{
				cBody = 0.5*(Foam::tanh(intSpan_*dist/cellSize) + 1.0);
			}
			else
			{
				cBody = 0.5*(-1.0*Foam::tanh(intSpan_*dist/cellSize) + 1.0);
			}
		}
	}

	// Clip field values
	lambda[cellI] += cBody;
	lambda[cellI] = min(max(0.0, lambda[cellI]), 1.0);

	return neighbours;
}

//---------------------------------------------------------------------------//

bool stlModel::isBodyInMesh()
{
	const boundBox ibBound(bounds());

	forAll(geometricD_, dir)
	{
		if (geometricD_[dir] == 1)
		{
			if
			(!(
				meshBounds_.max()[dir] >= ibBound.min()[dir]
			 && meshBounds_.min()[dir] <= ibBound.max()[dir]
			))
			{
				return false;
			}
		}
	}

	return true;
}

//---------------------------------------------------------------------------//

bool stlModel::isPointInBody
(
	point pointI
)
{
	pointField points(1, pointI);
	boolList returnList = triSurfSearch_().calcInside(points);

	return returnList[0];
}

//---------------------------------------------------------------------------//

label stlModel::findCellInBody()
{
	labelHashSet visited;

	const pointField& cellCenters = mesh_.C();

	if (cellToStart_ >= mesh_.nCells())
	{
		cellToStart_ = 0;
	}

	autoPtr<DynamicLabelList> pending(new DynamicLabelList(1, cellToStart_));
	autoPtr<DynamicLabelList> nextPending(new DynamicLabelList);

	label iterCount(0); const label iterMax(mesh_.nCells());

	while (pending().size() > 0 and iterCount < iterMax)
	{
		nextPending().clear();
		forAll(pending(), cellToCheck)
		{
			if (!visited.found(pending()[cellToCheck]))
			{
				visited.insert(pending()[cellToCheck]);
				iterCount++;

				if (isPointInBody(cellCenters[pending()[cellToCheck]]))
				{
					return pending()[cellToCheck];
				}
				else
				{
					nextPending().append
					(
						mesh_.cellCells()[pending()[cellToCheck]]
					);
				}
			}
		}

		autoPtr<DynamicLabelList> helperPtr(pending.ptr());
		pending.reset(nextPending.ptr());
		nextPending = std::move(helperPtr);
	}

	return -1;
}

void stlModel::findProcBoundaryCells
(
	label cellI,
	List<DynamicLabelList>& neighboursToSend
)
{
	forAll(mesh_.cells()[cellI], fI)
	{
		label faceI = mesh_.cells()[cellI][fI];

		if (!mesh_.isInternalFace(faceI))
		{
			label facePatchI
			(
				mesh_.boundaryMesh().whichPatch(faceI)
			);
			const polyPatch& patch
				= mesh_.boundaryMesh()[facePatchI];

			if (patch.type() ==	"processor")
			{
				const processorPolyPatch& procPatch
					= refCast<const processorPolyPatch>(patch);
				label iProc =
				(
					Pstream::myProcNo() == procPatch.myProcNo()
				)
					? procPatch.neighbProcNo()
					: procPatch.myProcNo();

				neighboursToSend[iProc].append(patch.whichFace(faceI));
			}
		}
	}
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void stlModel::generateLambda
(
	volScalarField& lambda
)
{
	internalCells_[Pstream::myProcNo()].clear();

	//- Check if body lies in mesh
	Info << "Checking if body intersects mesh" << endl;
	label pendingSize = 1;
	if (!isBodyInMesh())
	{
		pendingSize = 0;
	}

	// Octree traversal through mesh to find seed cell
	const pointField& cellCenters = mesh_.C();
	cellToStart_ = findCellInBody();

	if (cellToStart_ == -1)
	{
		pendingSize = 0;
	}

	// Octree traversal through mesh to determine lambda
	Field<label> visited(mesh_.nCells(), 0);
	autoPtr<DynamicLabelList> pending(new DynamicLabelList(pendingSize, cellToStart_));
	autoPtr<DynamicLabelList> nextPending(new DynamicLabelList);
	autoPtr<List<DynamicLabelList>> neighboursToSend
	(
		new List<DynamicLabelList>(Pstream::nProcs())
	);

	// Find the total number of empty patches
	label nEmptyDirs(0);
	forAll(mesh_.boundaryMesh(), patchI)
	{
		const polyPatch& patch = mesh_.boundaryMesh()[patchI];
		if (patch.type() == "empty")
		{
			nEmptyDirs++;
		}
	}

	HashTable<bool, label, Hash<label>> cellInside(128);

	label iterCount = 0; const label iterMax = mesh_.nCells();
	reduce(pendingSize, maxOp<label>());
	while (pendingSize > 0 && iterCount < iterMax)
	{
		nextPending().clear();

		forAll(pending(), cellToCheck)
		{
			const label cellI = pending()[cellToCheck];

			if (!cellInside.found(cellI))
			{
				iterCount++;

				if (isPointInBody(cellCenters[cellI]))
				{
					cellInside.set(cellI, true);

					const labelList& neighbours = mesh_.cellCells(cellI);
					nextPending().append(neighbours);

					label nProcFaces = mesh_.cells()[cellI].size();
					nProcFaces -= mesh_.cellCells()[cellI].size();
					nProcFaces -= nEmptyDirs;
					if (nProcFaces == 0)
					{
						continue;
					}

					findProcBoundaryCells(cellI, neighboursToSend());
				}
				else
				{
					cellInside.set(cellI, false);
				}
			}
		}

		// Send face indices to neighbours
		PstreamBuffers pBufsIFaces(Pstream::commsTypes::nonBlocking);
		for (label proci = 0; proci < Pstream::nProcs(); proci++)
		{
			if (proci != Pstream::myProcNo())
			{
				UOPstream sendIFaces(proci, pBufsIFaces);
				sendIFaces << neighboursToSend()[proci];
				neighboursToSend()[proci].clear();
			}
		}
		pBufsIFaces.finishedSends();

		// Recieve face indices and add to check
		for (label proci = 0; proci < Pstream::nProcs(); proci++)
		{
			if (proci != Pstream::myProcNo())
			{
				UIPstream recvIFaces(proci, pBufsIFaces);
				DynamicLabelList recvIFacesList(recvIFaces);

				// Find cells for faces
				forAll(recvIFacesList, rFace)
				{
					label faceI = recvIFacesList[rFace];

					// Find the cell
					forAll(mesh_.boundaryMesh(), patchI)
					{
						if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI]))
						{
							const processorPolyPatch& procPatch
								= refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

							// Get neighbouring processor id
							label iProc =
							(
								Pstream::myProcNo() == procPatch.myProcNo()
							)
								? procPatch.neighbProcNo()
								: procPatch.myProcNo();

							if (iProc == proci)
							{
								label rCellI = mesh_.boundaryMesh()[patchI].faceCells()[faceI];
								nextPending().append(rCellI);
							}
						}
					}
				}
			}
		}

		// Clear processor stream buffer
		pBufsIFaces.clear();

		// Clear pending queue and setup next wave
		autoPtr<DynamicLabelList> helperPtr(pending.ptr());
		pending.reset(nextPending.ptr());
		nextPending = std::move(helperPtr);

		// Check if all processors finished
		pendingSize = pending().size();
		reduce(pendingSize, maxOp<label>());
	}

	// Classify all cells found by octree
	forAll(cellInside.toc(), i)
	{
		label cellI = cellInside.toc()[i];
		bool centreInside = cellInside[cellI];

		pointField& points = cellPoints_[cellI];
		boolList verticesInside = triSurfSearch_().calcInside(points);

		classifyCell
		(
			cellI,
			centreInside,
			verticesInside,
			lambda
		);
	}

	// Update octree start cell for next call
	if (internalCells_[Pstream::myProcNo()].size() > 0)
	{
		cellToStart_ = min(internalCells_[Pstream::myProcNo()]);
	}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
