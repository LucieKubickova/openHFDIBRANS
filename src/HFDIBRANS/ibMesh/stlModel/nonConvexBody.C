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

#include "nonConvexBody.H"

#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "volFields.H"
#include "pointField.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonConvexBody::nonConvexBody
(
	const fvMesh& mesh,
	scalar thrSurf,
	scalar intSpan,
	bool sdBasedLambda,
	autoPtr<triSurfaceMesh>& triSurfMesh,
	autoPtr<triSurfaceSearch>& triSurfSearch
)
:
	stlModel
	(
		mesh,
		thrSurf,
		intSpan,
		sdBasedLambda,
		triSurfMesh,
		triSurfSearch
	)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

nonConvexBody::~nonConvexBody()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

labelList nonConvexBody::findBboxCells
(
	label cellToCheck,
	bool& isInsideBB,
	List<DynamicLabelList>& bBoxCells,
	Field<label>& visited
)
{
	labelList returnList;

	if (visited[cellToCheck] == 0)
	{
		visited[cellToCheck] = 1;
		vector cellCenter = mesh_.C()[cellToCheck];
		label partCheck(0);

		forAll(minBbox_, vecI)
		{
			if
			(
				cellCenter[vecI] >= minBbox_[vecI]
			 && cellCenter[vecI] <= maxBbox_[vecI]
			)
			{
				partCheck++;
			}
		}

		bool cellInside;
		if (partCheck == 3)
		{
			cellInside = true;
		}
		else
		{
			cellInside = false;
		}

		if (cellInside)
		{
			bBoxCells[Pstream::myProcNo()].append(cellToCheck);
			isInsideBB = true;
		}

		if (!isInsideBB || cellInside)
		{
			returnList = mesh_.cellCells()[cellToCheck];
		}
	}

	return returnList;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void nonConvexBody::generateLambda
(
	volScalarField& lambda
)
{
	internalCells_[Pstream::myProcNo()].clear();

	// Check if body lies in mesh
	Info<< "Checking if body intersects mesh" << endl;
	label pendingSize = 1;
	if (!isBodyInMesh())
	{
		pendingSize = 0;
	}

	// Reduce computational domain to the body bounding box
	scalar inflFact(2*sqrt(mesh_.magSf()[0]));
	vector unitVec(1, 1, 1);
	minBbox_ = bounds().min() - unitVec*inflFact;
	maxBbox_ = bounds().max() + unitVec*inflFact;

	// Octree traversal through mesh to determine lambda
	Field<label> visited(mesh_.nCells(), 0);
	autoPtr<DynamicLabelList> pending
	(
		new DynamicLabelList(pendingSize, cellToStart_)
	);
	autoPtr<DynamicLabelList> nextPending(new DynamicLabelList);
	List<DynamicLabelList> bBoxCells(Pstream::nProcs());
	autoPtr<List<DynamicLabelList>> neighboursToSend
	(
		new List<DynamicLabelList>(Pstream::nProcs())
	);

	// Find the total number of empty patches
	nEmptyDirs_ = 0;
	forAll(mesh_.boundaryMesh(), patchI)
	{
		const polyPatch& patch = mesh_.boundaryMesh()[patchI];
		if (isA<emptyPolyPatch>(patch))
		{
			nEmptyDirs_++;
		}
	}

	bool isInsideBB(false);
	label iterCount = 0; const label iterMax = mesh_.nCells();
	reduce(pendingSize, maxOp<label>());
	while (pendingSize > 0 && iterCount < iterMax)
	{
		iterCount++;
		nextPending().clear();

		forAll(pending(), cellToCheck)
		{
			const label cellI = pending()[cellToCheck];
			nextPending().append
			(
				findBboxCells
				(
					cellI,
					isInsideBB,
					bBoxCells,
					visited
				)
			);

			// Skip outside cells
			if (!isInsideBB)
			{
				continue;
			}

			// Skip if nNeighbours == nProcs
			label nProcFaces = mesh_.cells()[cellI].size();
			nProcFaces -= mesh_.cellCells()[cellI].size();
			nProcFaces -= nEmptyDirs_;
			if (nProcFaces == 0)
			{
				continue;
			}

			findProcBoundaryCells(cellI, neighboursToSend());
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
						const polyPatch& patch = mesh_.boundaryMesh()[patchI];
						if (isA<processorPolyPatch>(patch))
						{
							const processorPolyPatch& procPatch =
								refCast<const processorPolyPatch>(patch);

							// Get neighbouring processor id
							label iProc =
								(Pstream::myProcNo() == procPatch.myProcNo())
									? procPatch.neighbProcNo()
									: procPatch.myProcNo();

							if (iProc == proci)
							{
								label rCellI = patch.faceCells()[faceI];
								nextPending().append
								(
									findBboxCells
									(
										rCellI,
										isInsideBB,
										bBoxCells,
										visited
									)
								);
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

	// Find potent surface cells
	//DynamicLabelList potentSurfCells = findPotentSurfCells(lambda, cellInside);

	// Classify all cells found by octree
	Info<< "Calculating lambda values" << endl;
	const vector sdSpan(4.0*(mesh_.bounds().max() - mesh_.bounds().min()));
	forAll(bBoxCells[Pstream::myProcNo()], i)
	{
		label cellI = bBoxCells[Pstream::myProcNo()][i];
		bool centerInside = isPointInBody(mesh_.C()[cellI]);

		classifyCell(lambda, sdSpan, cellI, centerInside);
	}

	// Update octree start cell for next call
	// if (internalCells_[Pstream::myProcNo()].size() > 0)
	// {
	// 	cellToStart_ = min(internalCells_[Pstream::myProcNo()]);
	// }
}


// ************************************************************************* //
