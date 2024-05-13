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

InNamspace
    Foam

Description
    implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with RAS turbulence modeling approach and
    wall functions (Kubickova and Isoz, 2023)

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef HFDIBKOmegaSST_H
#define HFDIBKOmegaSST_H

#include "HFDIBRASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBRASModels
{

/*---------------------------------------------------------------------------*\
                           Class HFDIBKOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicHFDIBTurbulenceModel>
class HFDIBKOmegaSST
:
    public eddyViscosity<HFDIBRASModel<BasicHFDIBTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        HFDIBKOmegaSST(const HFDIBKOmegaSST&) = delete;

        //- No copy assignment
        void operator=(const HFDIBKOmegaSST&) = delete;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar b1_;
            dimensionedScalar c1_;

            //- Flag to include the F3 term
            Switch F3_;


        // Fields

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

            //- Turbulent kinetic energy field [m^2/s^2]
            volScalarField k_;

            //- Specific dissipation rate field [1/s]
            volScalarField omega_;

        // Decay control

            //- Flag to include the decay control
            Switch decayControl_;
            dimensionedScalar kInf_;
            dimensionedScalar omegaInf_;

        // HFDIBRANS
	        volScalarField lambda_;
            volScalarField kSurface_;
            volScalarField omegaGSurface_;
	        volScalarField ki_;
            volScalarField kQ_;
	        volScalarField nu_;

        // Dictionaries

            IOdictionary HFDIBDEMDict_;

        // Equations settings
            word kSurfaceType_;
            word omegaGSurfaceType_;
            scalar kBoundaryValue_;
            scalar omegaGBoundaryValue_;
            scalar tolKEqn_;
            label maxKEqnIters_;

        // Options
            bool useKQ_;

    // Protected Member Functions

        //- Set decay control with kInf and omegaInf
        void setDecayControl(const dictionary& dict);

        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        virtual tmp<volScalarField> F2() const;
        virtual tmp<volScalarField> F3() const;
        virtual tmp<volScalarField> F23() const;

        //- Return the blended field
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        //- Return the internal blended field
        tmp<volScalarField::Internal> blend
        (
            const volScalarField::Internal& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK(const volScalarField& F1) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega(const volScalarField& F1) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField::Internal> beta
        (
            const volScalarField::Internal& F1
        ) const
        {
            return tmp<volScalarField::Internal>::New
            (
                this->type() + ":beta",
                blend(F1, beta1_, beta2_)
            );
        }

        tmp<volScalarField::Internal> gamma
        (
            const volScalarField::Internal& F1
        ) const
        {
            return tmp<volScalarField::Internal>::New
            (
                this->type() + ":gamma",
                blend(F1, gamma1_, gamma2_)
            );
        }

        virtual void correctNut(const volScalarField& S2);

        virtual void correctNut();

        //- Return square of strain rate
        virtual tmp<volScalarField> S2
        (
            const volTensorField& gradU
        ) const;

        //- Return k production rate
        virtual tmp<volScalarField::Internal> Pk
        (
            const volScalarField::Internal& G
        ) const;

        //- Return epsilon/k which for standard RAS is betaStar*omega
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU
        ) const;

        //- Return (G/nu)_0
        virtual tmp<volScalarField::Internal> GbyNu0
        (
            const volTensorField& gradU,
            const volScalarField& S2
        ) const;

        //- Return G/nu
        virtual tmp<volScalarField::Internal> GbyNu
        (
            const volScalarField::Internal& GbyNu0,
            const volScalarField::Internal& F2,
            const volScalarField::Internal& S2
        ) const;

        virtual tmp<fvScalarMatrix> kSource() const;

        virtual tmp<fvScalarMatrix> omegaSource() const;

        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField::Internal& S2,
            const volScalarField::Internal& gamma,
            const volScalarField::Internal& beta
        ) const;

        virtual void matrixManipulate(fvScalarMatrix& eqn, volScalarField& phi, volScalarField& surface);

public:

    typedef typename BasicHFDIBTurbulenceModel::alphaField alphaField;
    typedef typename BasicHFDIBTurbulenceModel::rhoField rhoField;
    typedef typename BasicHFDIBTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("HFDIBKOmegaSST");


    // Constructors

        //- Construct from components
        HFDIBKOmegaSST
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = HFDIBturbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~HFDIBKOmegaSST() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>::New
            (
                "DkEff",
                alphaK(F1)*this->nut_ + this->nu()
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>::New
            (
                "DomegaEff",
                alphaOmega(F1)*this->nut_ + this->nu()
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- second image of the correct function to accept arguments
        virtual void correct(openHFDIBRANS& HFDIBRANS);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBRASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBKOmegaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
