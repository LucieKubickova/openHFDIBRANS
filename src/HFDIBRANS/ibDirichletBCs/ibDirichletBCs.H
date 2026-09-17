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
    Foam::ibDirichletBCs

Description
    Implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with RAS turbulence modeling approach and
    wall functions (Kubickova and Isoz, 2023)

SourceFiles
    ibDirichletBCs.C
    ibDirichletBCsI.H

Contributors
	Federico Municchi (2016),
	Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),

\*---------------------------------------------------------------------------*/

#ifndef ibDirichletBCs_H
#define ibDirichletBCs_H

#include "nearWallDist.H"

#include "ibMesh.H"
#include "boundaryCell.H"
#include "surfaceCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibDirichletBCs Declaration
\*---------------------------------------------------------------------------*/

class ibDirichletBCs
{

private:

    // Private Data

		//- fvMesh
		const fvMesh& mesh_;

		//- ibMesh
		ibMesh& ibMesh_;

		//- Lambda field
		const volScalarField& body_;

		//- Options
		scalar thrSurf_;
		bool useYEff_;
		word uTauType_;
		scalar uTauCoeff_;

		//- Fields
		List<List<scalar>> nutAtIB_;
		List<List<scalar>> kAtIB_;
		List<List<scalar>> uTauAtIB_;

		//- References to shared data
		List<DynamicList<boundaryCell>>& boundaryCells_;
		List<DynamicList<surfaceCell>>& surfaceCells_;
		labelField& isBoundaryCell_;

		//- Dictionaries
		IOdictionary turbulenceProperties_;
		IOdictionary HFDIBDEMDict_;
		IOdictionary fvSchemes_;
		dictionary HFDIBBCsDict_;

		//- Simulation type
		word simulationType_;

		//- Saves
		volScalarField yPlusi_;
		volScalarField uTaui_;
		volScalarField nuti_;

		//- Boundary conditions
		word UBC_;

		//- RAS wall function names
		word nutWF_;
		word kWF_;
		word omegaWF_;
		word epsilonWF_;

		//- RAS parameters
		const scalar kappa_;
		const scalar E_;
		const scalar Cmu_;
		const scalar Ceps2_;
		const scalar beta1_;
		scalar Cmu75_;
		scalar Cmu25_;
		scalar Cmu5_;
		scalar yPlusLam_;

    // Private Member Functions

		//- Calculate laminar y+
		void calcYPlusLam();

		//- Get zero scalar
		inline scalar ibZero(scalar dummy);

		//- Get zero vector
		inline vector ibZero(vector dummy);

		//- Get unit scalar
		inline scalar ibOne(scalar dummy);

		//- Get unit vector
		inline vector ibOne(vector dummy);


public:

    // Constructors

		//- Contruct from components
		ibDirichletBCs
		(
			const fvMesh& mesh,
			ibMesh& ibMesh,
			const volScalarField& body,
			List<DynamicList<boundaryCell>>& boundaryCells,
			List<DynamicList<surfaceCell>>& surfaceCells,
			labelField& isBoundaryCell
		);

    //- Destructor
    ~ibDirichletBCs();

    // Member Functions

		//- Set size to list
		void setSizeToLists();

		//- Velocity boundary condition
		void UAtIB(List<vector>& UIB, volVectorField& U);

		//- Scalar boundary condition
		void TAtIB(List<scalar>& TIB, scalar TIn);

		//- Calculate friction velocity
		void updateUTauAtIB(volScalarField& k);

		//- Correct nut at the immersed boundary
		void nutAtIB(volScalarField& k, volScalarField& nu);

		//- k boundary condition
		void kAtIB
		(
			List<scalar>& kIB,
			volScalarField& k,
			volScalarField& nu
		);

		//- omegaG boundary condition
		void omegaGAtIB
		(
			List<scalar>& omegaIB,
			List<scalar>& GIB,
			volScalarField::Internal& G,
			const volVectorField& U,
			volScalarField& k,
			volScalarField& nu
		);

		//- epsilonG boundary condition
		void epsilonGAtIB
		(
			List<scalar>& epsilonIB,
			List<scalar>& GIB,
			volScalarField::Internal& G,
			const volVectorField& U,
			volScalarField& k,
			volScalarField& nu
		);

		//- Calculate wall shear stress
		void calculateWallShearStress(volVectorField& tauw, const volVectorField& U, const volScalarField& nu);

		//- Calculate forces
		void calculateForces(volVectorField& fN, volVectorField& fT, volVectorField& tauw, const volScalarField& p, dictionary forceDict);

		//- Calculate force coefficients
		void calculateForceCoeffs(scalar& Cl, scalar& Cd, volVectorField& fN, volVectorField& fT, dictionary forceDict);

		//- Save friction velocity
		void saveUTau();

		//- Calcualte surface gradient
		void snGradUAtIB(const volVectorField& U, List<DynamicList<vector>>& snGradU);

		// Access

			//- Get immersed yPlus field
			inline volScalarField& getYPlusi();

			//- Get yPlusLam
			inline scalar getYPlusLam();

			//- Get E coefficient
			inline scalar getE();

			//- Get nut at immersed boundary
			inline List<List<scalar>>& getNutAtIB();

			//- Get k at immersed boundary
			inline List<List<scalar>>& getKAtIB();

			//- Get imposed nut field
			inline volScalarField& getNuti();

			//- Get imposed friction velocity field
			inline volScalarField& getUTaui();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ibDirichletBCsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
