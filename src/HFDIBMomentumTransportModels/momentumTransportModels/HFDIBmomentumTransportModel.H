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

#ifndef HFDIBmomentumTransportModel_H
#define HFDIBmomentumTransportModel_H

#include "IOdictionary.H"
#include "primitiveFieldsFwd.H"
#include "volFields.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "geometricOneField.H"
#include "momentumTransportModel.H"
#include "viscosity.H"
#include "openHFDIBRANS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class HFDIBmomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class HFDIBmomentumTransportModel
:
    //~ public IOdictionary
    public momentumTransportModel
{
protected:

    // Protected data

        //~ const Time& runTime_;
        //~ const fvMesh& mesh_;

        //~ const volVectorField& U_;
        //~ const surfaceScalarField& alphaRhoPhi_;
        //~ const surfaceScalarField& phi_;

        //~ const viscosity& viscosity_;


    // Protected member functions

        //~ static typeIOobject<IOdictionary> readModelDict
        //~ (
        //~     const objectRegistry& obr,
        //~     const word& group,
        //~     bool registerObject = false
        //~ );


        //~ template<class HFDIBMomentumTransportModel>
        //~ static inline autoPtr<HFDIBMomentumTransportModel> New
        //~ (
            //~ const typename HFDIBMomentumTransportModel::alphaField& alpha,
            //~ const typename HFDIBMomentumTransportModel::rhoField& rho,
            //~ const volVectorField& U,
            //~ const surfaceScalarField& alphaRhoPhi,
            //~ const surfaceScalarField& phi,
            //~ const viscosity& viscosity
        //~ );


public:

    //- Runtime type information
    TypeName("HFDIBmomentumTransport");


    // Constructors

        //- Construct from components
        HFDIBmomentumTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );

        //- Disallow default bitwise copy construction
        HFDIBmomentumTransportModel(const HFDIBmomentumTransportModel&) = delete;


    //- Destructor
    virtual ~HFDIBmomentumTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read() = 0;

        const Time& time() const
        {
            return runTime_;
        }

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const = 0;

        //- Helper function to return the name of the turbulence G field
        inline word GName() const
        {
            return modelName("G");
        }

        //- Access function to velocity field
        inline const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to phase flux field
        inline const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Return the volumetric flux field
        virtual tmp<surfaceScalarField> phi() const;

        //- Access function to fluid properties
        const class Foam::viscosity& properties() const
        {
            return viscosity_;
        }

        //- Return the near wall distance
        const volScalarField::Boundary& y() const;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return this->viscosity_.nu();
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return this->viscosity_.nu(patchi);
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const = 0;

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const = 0;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const = 0;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const = 0;

        //- Return the stress tensor [m^2/s^2]
        virtual tmp<volSymmTensorField> sigma() const = 0;

        //- Validate the turbulence fields after construction
        //  Update derived fields as required
        virtual void validate();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct(openHFDIBRANS& HFDIBRANS) = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HFDIBmomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBmomentumTransportModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
