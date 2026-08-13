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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),
	Vít Večerník (2026-*)

\*---------------------------------------------------------------------------*/

#include "HashTable.H"
#include "Pstream.H"
#include "labelList.H"
#include "pointIndexHit.H"
#include "vector.H"
#include <algorithm>

#include "stlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stlModel::stlModel
(
	const fvMesh& mesh,
	scalar& thrSurf,
	scalar& intSpan,
	autoPtr<triSurfaceMesh>& triSurfMesh,
	autoPtr<triSurfaceSearch>& triSurfSearch
)
:
	mesh_(mesh),
	thrSurf_(thrSurf),
	intSpan_(intSpan),
	triSurfMesh_(triSurfMesh),
	triSurfSearch_(triSurfSearch),
	geometricD_(mesh_.geometricD()),
	meshBounds_(mesh_.points(), false),
	cellToStart_(0),
	cellPoints_(mesh_.nCells()),
	lastCellPoints_(0)
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
	bool& insideIB,
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
		else
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
				Info << "stlModel: missed nearest point for cell "
					 << cellI << endl;
			}

			const scalar cellSize = Foam::pow(mesh_.V()[cellI], 0.333);

			if (centreInside)
			{
				cBody = 0.5*(Foam::tanh(intSpan_*dist/cellSize) + 1.0);
			}
			else
			{
				cBody = 0.5*(-Foam::tanh(intSpan_*dist/cellSize) + 1.0);
			}
		}

		insideIB = true;
		neighbours = mesh_.cellCells()[cellI];
	}
	else if (!insideIB)
	{
		neighbours = mesh_.cellCells()[cellI];
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
	bool insideMesh = false;
	const pointField ibCorners(ibBound.points());

	forAll(ibCorners, pI)
	{
		bool cornerOk = true;
		forAll(geometricD_, dir)
		{
			if (geometricD_[dir] == 1)
			{
				if
				(
					ibCorners[pI][dir] <= meshBounds_.min()[dir]
				 || ibCorners[pI][dir] >= meshBounds_.max()[dir]
				)
				{
					cornerOk = false;
				}
			}
		}
		if (cornerOk)
		{
			insideMesh = true;
			break;
		}
	}

	return insideMesh;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void stlModel::generateLambda
(
	volScalarField& lambda
)
{
	internalCells_[Pstream::myProcNo()].clear();

	//- Check if atleast one body corner lies in mesh
	bool insideMesh = isBodyInMesh();

	if (!insideMesh)
	{
		Info << "stlModel: body bounding box lies outside the mesh. "
			 << "Aborting lambda generation." << endl;
		return;
	}

	// Octree traversal through mesh cells
	if (cellToStart_ >= mesh_.nCells())
	{
		cellToStart_ = 0;
	}

	Field<label> visited(mesh_.nCells(), 0);
	labelList pending(1, cellToStart_);
	label iterCount = 0;
	const label iterMax = mesh_.nCells();

	const boundBox ibBound(bounds());
	bool insideIB = false;
	bool insideIBBound = false;
	HashTable<pointField, label, Hash<label>> newCellPtsCache(0);
	DynamicLabelList nextPending;

	while (pending.size() > 0 && iterCount < iterMax)
	{
		iterCount++;
		nextPending.clear();

		forAll(pending, pi)
		{
			const label cellI = pending[pi];

			if (visited[cellI]) continue;
			visited[cellI] = 1;

			// Get vertex positions, use cache whenever possible
			pointField pts;
			if (lastCellPoints_.found(cellI))
			{
				pts = lastCellPoints_[cellI];
			}
			else
			{
				pts = cellPoints_[cellI];
			}

			// Check if any vertex lies inside the IB box
			bool inBB = false;
			forAll(pts, pI)
			{
				if (ibBound.contains(pts[pI]))
				{
					inBB = true;
					break;
				}
			}

			if (inBB)
			{
				insideIBBound = true;
				cellToStart_ = cellI;
				newCellPtsCache.insert(cellI, pts);

				boolList verticesInside = triSurfSearch_().calcInside(pts);
                const pointField centreField(1, mesh_.C()[cellI]);
				bool centreInside = triSurfSearch_().calcInside(centreField)[0];

				if
				(
					std::any_of
					(
						verticesInside.begin(), verticesInside.end(),
						[](bool b){ return b; }
					)
				 || centreInside
				 || !insideIB
				)
				{
					const labelList neighbours
					(
						classifyCell
						(
							cellI,
							insideIB,
							centreInside,
							verticesInside,
							lambda
						)
					);
					nextPending.append(neighbours);
				}
			}
			else if (!insideIB && !insideIBBound)
			{
				nextPending.append(mesh_.cellCells()[cellI]);
			}
		}
		pending = nextPending;
	}

	// Update cached cell-point positions
	lastCellPoints_ = newCellPtsCache;

	// Update octree start cell for next call
	if (internalCells_[Pstream::myProcNo()].size() > 0)
	{
		cellToStart_ = min(internalCells_[Pstream::myProcNo()]);
	}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
