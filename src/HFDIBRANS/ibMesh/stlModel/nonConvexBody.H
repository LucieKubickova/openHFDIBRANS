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
	Foam::nonConvexBody

Description
	Class for the generation of a lambda field based on an STL file.

SourceFiles
	nonConvexBody.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),
	Vít Večerník (2026-*)

\*---------------------------------------------------------------------------*/

#ifndef nonConvexBody_H
#define nonConvexBody_H

#include "stlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class nonConvexBody Declaration
\*---------------------------------------------------------------------------*/

class nonConvexBody
:
	public stlModel
{

private:

	typedef DynamicList<label> DynamicLabelList;

	// Private Data

		//- Bounding box
		vector minBbox_;
		vector maxBbox_;

	// Private Member Functions

		// //- Find potent surface cells
		// DynamicLabelList findPotentSurfCells
		// (
		// 	volScalarField& lambda,
		// 	HashTable<bool, label, Hash<label>>& cellInside
		// );

		//- Find cells inside boundary box
		labelList findBboxCells
		(
			label cellToCheck,
			bool& isInsideBB,
			List<DynamicLabelList>& bBoxCells,
			Field<label>& visited
		);


public:

	// Constructors

		nonConvexBody
		(
			const fvMesh& mesh,
			scalar thrSurf,
			scalar intSpan,
			bool sdBasedLambda_,
			autoPtr<triSurfaceMesh>& triSurfMesh,
			autoPtr<triSurfaceSearch>& triSurfSearch
		);


	//- Destructor
	~nonConvexBody();


	// Member Functions

		//- Generate the lambda field based on an STL
		void generateLambda(volScalarField& lambda) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
