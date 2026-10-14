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

#include "IncompressibleHFDIBTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::
IncompressibleHFDIBTurbulenceModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
:
    HFDIBTurbulenceModel
    <
        geometricOneField,
        geometricOneField,
        incompressibleHFDIBTurbulenceModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<IncompressibleHFDIBTurbulenceModel>
    (
        static_cast<IncompressibleHFDIBTurbulenceModel*>(
        HFDIBTurbulenceModel
        <
            geometricOneField,
            geometricOneField,
            incompressibleHFDIBTurbulenceModel,
            TransportModel
        >::New
        (
            geometricOneField(),
            geometricOneField(),
            U,
            phi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::devReff
(
    const volVectorField& U
) const
{
    return devRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::devRhoReff() const
{
    NotImplemented;

    return devReff();
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::devRhoReff
(
    const volVectorField& U
) const
{
    NotImplemented;

    return nullptr;
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleHFDIBTurbulenceModel<TransportModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}


// ************************************************************************* //
