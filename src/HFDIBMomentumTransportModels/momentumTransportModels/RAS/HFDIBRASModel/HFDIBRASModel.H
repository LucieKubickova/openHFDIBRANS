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

#ifndef HFDIBRASModel_H
#define HFDIBRASModel_H

#include "HFDIBmomentumTransportModel.H"
#include "generalisedNewtonianViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HFDIBRASModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicHFDIBMomentumTransportModel>
class HFDIBRASModel
:
    public BasicHFDIBMomentumTransportModel
{

protected:

    // Protected data

        //- HFDIBRAS coefficients dictionary
        dictionary HFDIBRASDict_;

        //- Turbulence on/off flag
        Switch turbulence_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Lower limit of k
        dimensionedScalar kMin_;

        //- Lower limit of epsilon
        dimensionedScalar epsilonMin_;

        //- Lower limit for omega
        dimensionedScalar omegaMin_;

        //- Run-time selectable generalised Newtonian viscosity model
        autoPtr<laminarModels::generalisedNewtonianViscosityModel>
            viscosityModel_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs(const word& type);


public:

    typedef typename BasicHFDIBMomentumTransportModel::alphaField alphaField;
    typedef typename BasicHFDIBMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("HFDIBRAS");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            HFDIBRASModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const viscosity& viscosity
            ),
            (alpha, rho, U, alphaRhoPhi, phi, viscosity)
        );


    // Constructors

        //- Construct from components
        HFDIBRASModel
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );

        //- Disallow default bitwise copy construction
        HFDIBRASModel(const HFDIBRASModel&) = delete;


    // Selectors

        //- Return a reference to the selected HFDIBRAS model
        static autoPtr<HFDIBRASModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~HFDIBRASModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the lower allowable limit for k (default: small)
        const dimensionedScalar& kMin() const
        {
            return kMin_;
        }

        //- Return the lower allowable limit for epsilon (default: small)
        const dimensionedScalar& epsilonMin() const
        {
            return epsilonMin_;
        }

        //- Return the lower allowable limit for omega (default: small)
        const dimensionedScalar& omegaMin() const
        {
            return omegaMin_;
        }

        //- Allow kMin to be changed
        dimensionedScalar& kMin()
        {
            return kMin_;
        }

        //- Allow epsilonMin to be changed
        dimensionedScalar& epsilonMin()
        {
            return epsilonMin_;
        }

        //- Allow omegaMin to be changed
        dimensionedScalar& omegaMin()
        {
            return omegaMin_;
        }

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return viscosityModel_->nu();
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return viscosityModel_->nu(patchi);
        }

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return volScalarField::New
            (
                IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
                this->nut() + this->nu()
            );
        }

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return this->nut(patchi) + this->nu(patchi);
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- second image of the correct function to accept arguments
        virtual void correct(openHFDIBRANS& HFDIBRANS);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HFDIBRASModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBRASModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
