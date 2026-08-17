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

		//- Minimum contribution to treat as inside
		const scalar& thrSurf_;

		//- Width parameter for signed distance lambda
		const scalar& intSpan_;

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
		List<pointField> cellPoints_;

		//- Cache of cell-vertex positions
		HashTable<pointField, label, Hash<label>> lastCellPoints_;

		//- List of per processor internal cells
		List<DynamicLabelList> internalCells_;


	// Private Member Functions

		//- Process cell, clasify it and return neighbors
		labelList classifyCell
		(
			label cellI,
			bool& insideIB,
			bool& centerInside,
			boolList& verticesInside,
			volScalarField& lambda
		);

		//- Check if STL body lies inside mesh
		bool isBodyInMesh();

public:

	// Constructors

		stlModel
		(
			const fvMesh& mesh,
			scalar& thrSurf,
			scalar& intSpan,
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

#endif

// ************************************************************************* //
