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

#ifndef incompressibleHFDIBMomentumTransportModel_H
#define incompressibleHFDIBMomentumTransportModel_H

#include "momentumTransportModel.H"
#include "HFDIBmomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class incompressibleHFDIBMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class incompressibleHFDIBMomentumTransportModel
:
    public HFDIBmomentumTransportModel
{

public:

    typedef geometricOneField alphaField;
    typedef geometricOneField rhoField;


protected:

    // Protected data

        alphaField alpha_;
        rhoField rho_;


public:

    //- Runtime type information
    // TypeName(HFDIBmomentumTransportModel::typeName_());


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            incompressibleHFDIBMomentumTransportModel,
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

        //- Construct
        incompressibleHFDIBMomentumTransportModel
        (
            const word& type,
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<incompressibleHFDIBMomentumTransportModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~incompressibleHFDIBMomentumTransportModel()
    {}


    // Member Functions

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return alpha_;
        }

        //- Return the density field
        const rhoField& rho() const
        {
            return rho_;
        }

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devSigma() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevSigma(volVectorField& U) const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devTau() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;
};


namespace incompressible
{
    typedef incompressibleHFDIBMomentumTransportModel HFDIBmomentumTransportModel;

    template<class IncompressibleHFDIBMomentumTransportModel>
    autoPtr<IncompressibleHFDIBMomentumTransportModel> New
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const viscosity& viscosity
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "incompressibleHFDIBMomentumTransportModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
