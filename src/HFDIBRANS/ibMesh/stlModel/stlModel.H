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

Class
	Foam::stlModel

Description
	Class for the generation of a lambda field based on an STL file.

SourceFiles
	stlModel.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),
	Vít Večerník (2026-*)

\*---------------------------------------------------------------------------*/

#ifndef stlModel_H
#define stlModel_H

#include "DynamicList.H"
#include "autoPtr.H"
#include "boolList.H"
#include "boundBox.H"
#include "fvMesh.H"
#include "pointField.H"
#include "triSurfaceMesh.H"
#include "triSurfaceSearch.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class stlModel Declaration
\*---------------------------------------------------------------------------*/

class stlModel
{

private:

	typedef DynamicList<label> DynamicLabelList;

	// Private Data

		//- Mesh
		const fvMesh& mesh_;

		//- Options
		const scalar thrSurf_;
		const scalar intSpan_;
		const bool sdBasedLambda_;

		//- STL surface mesh
		const autoPtr<triSurfaceMesh>& triSurfMesh_;
		const autoPtr<triSurfaceSearch>& triSurfSearch_;

		//- Active geometric directions
		const Vector<label> geometricD_;

		//- Mesh bounding box
		const boundBox meshBounds_;

		//- BFS starting cell
		label cellToStart_;

		//- List of all cell-certex positions
		List<labelList> cellPoints_;

		//- List of per processor internal cells
		List<DynamicLabelList> internalCells_;

		//- Number of empty directions
		label nEmptyDirs_;


	// Private Member Functions

		//- Find potent surface cells
		DynamicLabelList findPotentSurfCells
		(
			volScalarField& lambda,
			HashTable<bool, label, Hash<label>>& cellInside
		);

		//- Calculate lambda value for candidate cell
		void classifyCell
		(
			volScalarField& lambda,
			vector sdSpan,
			label cellI,
			bool centerInside
		);

		//- Check if STL body lies inside mesh
		bool isBodyInMesh();

		//- Check if point lies inside STL body
		bool isPointInBody(point pointI);

		//- Find a cell inside body
		label findCellInBody();

		//- Find processor boundary cells
		void findProcBoundaryCells
		(
			label cellI,
			List<DynamicLabelList>& neighborsToSend
		);


public:

	// Constructors

		stlModel
		(
			const fvMesh& mesh,
			scalar thrSurf,
			scalar intSpan,
			bool sdBasedLambda_,
			autoPtr<triSurfaceMesh>& triSurfMesh,
			autoPtr<triSurfaceSearch>& triSurfSearch
		);


	//- Destructor
	~stlModel();


	// Member Functions

		//- Generate the lambda field based on an STL
		void generateLambda(volScalarField& lambda);

		//- Return bounding box of the STL
		boundBox bounds() const
		{
			return boundBox(triSurfMesh_().points());
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
