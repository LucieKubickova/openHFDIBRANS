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
	Foam::openHFDIBRANS

Description
    Implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with RAS turbulence modeling approach and
    wall functions (Kubickova and Isoz, 2023).

SourceFiles
	openHFDIBRANS.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)

\*---------------------------------------------------------------------------*/

#ifndef openHFDIBRANS_H
#define openHFDIBRANS_H

#include "ibInterpolation.H"
#include "ibDirichletBCs.H"

#define ORDER 2

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class openHFDIBRANS Declaration
\*---------------------------------------------------------------------------*/

class openHFDIBRANS
{

private:

    // Private Data

		//- fvMesh
		const fvMesh& mesh_;

		//- Lambda field
		volScalarField& body_;

		//- Threshold surface
		scalar thrSurf_;

		//- Dictionaries
		IOdictionary HFDIBDEMDict_;
		IOdictionary fvSchemes_;
		dictionary HFDIBOuterSchemes_;

		//- Lists of cells relevant to imposed fields
		List<DynamicList<boundaryCell>> boundaryCells_;
		List<DynamicList<surfaceCell>> surfaceCells_;
		List<DynamicList<label>> internalCells_;
		labelField isBoundaryCell_;

		//- Mesh class with auxiliar functions
		autoPtr<ibMesh> ibMesh_;

		//- Interpolation class
		autoPtr<ibInterpolation> ibInterpolation_;

		//- Wall functions class
		autoPtr<ibDirichletBCs> ibDirichletBCs_;

		//- Options
		bool save_;
		bool cpDisToInner_;
		bool scaleDisG_;
		bool scaleG_;
		scalar scaleCoeff_;
		scalar scaleCoeffG_;
		bool useYEff_;
		bool assignNut_;

		//- Save
		word outDir_;


public:

    // Constructors

		//- Construct from components
		openHFDIBRANS(const fvMesh& mesh, volScalarField& body);

    //- Destructor
    ~openHFDIBRANS();

    // Member Functions

		//- Compute imposed U field
		void computeUi
		(
			volVectorField& U,
			volVectorField& Ui,
			word surfType
		);

		//- Compute imposed k field
		void computeKi
		(
			volScalarField& k,
			volScalarField& ki,
			volScalarField& nu
		);

		//- Compute imposed T field
		void computeTi
		(
			volScalarField& T,
			volScalarField& Ti,
			word surfType,
			scalar TIn
		);

		//- Update friction velocity
		void updateUTau(volScalarField& k);

		//- Correct nut values
		void correctNut
		(
			volScalarField& nut,
			volScalarField& k,
			volScalarField& nu
		);

		//- Correct omegaG values
		void correctOmegaG
		(
			volScalarField& omega,
			volScalarField::Internal& G,
			const volVectorField& U,
			volScalarField& k,
			volScalarField& nu,
			volScalarField& surface
		);

		//- Correct epsilonG values
		void correctEpsilonG
		(
			volScalarField& epsilon,
			volScalarField::Internal& G,
			const volVectorField& U,
			volScalarField& k,
			volScalarField& nu,
			volScalarField& surface
		);

		//- Create and return base surface field
		void createBaseSurface
		(
			volScalarField& surface,
			word surfType,
			scalar boundaryValue
		);

		//- Update surface field
		void updateSurface(volScalarField& surface, word surfType);

		//- Correct surface field by normal
		void correctSurfaceByNormal
		(
			volSymmTensorField& normSurface,
			volScalarField& surface,
			scalar bodyOnLimit
		);

		//- Calculate surface area
		void calculateSurfaceArea();

		//- Calculate wall shear stress
		void calculateWallShearStress
		(
			volVectorField& tauw,
			const volVectorField& U,
			const volScalarField& nu
		);

		//- Calculate forces
		void calculateForces
		(
			volVectorField& fN,
			volVectorField& fT,
			volVectorField& tauw,
			const volScalarField& p,
			dictionary forceDict
		);

		//- Calculate force coefficients
		void calculateForceCoeffs
		(
			scalar& Cl,
			scalar& Cd,
			volVectorField& fN,
			volVectorField& fT,
			dictionary forceDict
		);

		//- Correct y for blending functions in kOmegaSST
		void correctY(volScalarField& y, bool recreate);

		//- Bound fields
		void bound(volScalarField& phi, dimensionedScalar& phiMin);

		// Access

			//- Get imposed turbulent viscosity
			volScalarField& getNuti()
			{
				return ibDirichletBCs_->getNuti();
			}

			//- Get imposed friction velocity
			volScalarField& getUTaui()
			{
				return ibDirichletBCs_->getUTaui();
			}

			//- Get lambda field
			volScalarField& getLambda() const
			{
				return body_;
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
