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
	Foam::convexBody

Description
	Class for the generation of a lambda field based on an STL file of a convex
	body.

SourceFiles
	convexBody.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),
	Vít Večerník (2026-*)

\*---------------------------------------------------------------------------*/

#ifndef convexBody_H
#define convexBody_H

#include "stlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class convexBody Declaration
\*---------------------------------------------------------------------------*/

class convexBody
:
	public stlModel
{

private:

	typedef DynamicList<label> DynamicLabelList;

	// Private Member Functions

		//- Find a cell inside body
		label findCellInBody();

		//- Find potent surface cells
		DynamicLabelList findPotentSurfCells
		(
			volScalarField& lambda,
			HashTable<bool, label, Hash<label>>& cellInside
		);


public:

	// Constructors

		//- Construct from components
		convexBody
		(
			const fvMesh& mesh,
			scalar thrSurf,
			scalar intSpan,
			bool sdBasedLambda,
			autoPtr<triSurfaceMesh>& triSurfMesh,
			autoPtr<triSurfaceSearch>& triSurfSearch
		);


	//- Destructor
	virtual ~convexBody();


	// Member Functions

		//- Generate the lambda field based on an STL
		virtual void generateLambda(volScalarField& lambda) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
