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

#include "stlModel.H"

#include "volFields.H"
#include "pointFields.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

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

	cellPoints_.setSize(mesh_.nCells());
	forAll(mesh.C(), cellI)
	{
		cellPoints_[cellI] = mesh_.cellPoints()[cellI];
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

stlModel::~stlModel()
{}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void stlModel::classifyCell
(
	volScalarField& lambda,
	vector sdSpan,
	label cellI,
	bool centerInside
)
{
	const labelList& cellVertices = cellPoints_[cellI];
	// Vertex weight contribution:
	// each vertex = 0.5/N, centre = 0.5 -> 1.0 in total
	const scalar vertexWeight = 0.5/cellVertices.size();
	scalar cBody = 0;

	forAll(cellVertices, vertI)
	{
		bool vertexInside = isPointInBody(mesh_.points()[cellVertices[vertI]]);
		if (vertexInside)
		{
			cBody += vertexWeight;
		}
	}
	if (centerInside)
	{
		cBody += 0.5;
	}

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
			pointIndexHit hit
			(
				triSurfSearch_().nearest(mesh_.C()[cellI], sdSpan)
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

			if (centerInside)
			{
				cBody = 0.5*(Foam::tanh(intSpan_*dist/cellSize) + 1.0);
			}
			else
			{
				cBody = 0.5*(-1.0*Foam::tanh(intSpan_*dist/cellSize) + 1.0);
			}
		}

		// Clip field values
		lambda[cellI] += cBody;
		lambda[cellI] = min(max(0.0, lambda[cellI]), 1.0);
	}
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
			label facePatchI = mesh_.boundaryMesh().whichPatch(faceI);
			const polyPatch& patch = mesh_.boundaryMesh()[facePatchI];

			if (isA<processorPolyPatch>(patch))
			{
				const processorPolyPatch& procPatch
					= refCast<const processorPolyPatch>(patch);
				label iProc = (Pstream::myProcNo() == procPatch.myProcNo())
					? procPatch.neighbProcNo() : procPatch.myProcNo();

				neighboursToSend[iProc].append(patch.whichFace(faceI));
			}
		}
	}
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

boundBox stlModel::bounds() const
{
	return boundBox(triSurfMesh_().points());
}


// ************************************************************************* //
