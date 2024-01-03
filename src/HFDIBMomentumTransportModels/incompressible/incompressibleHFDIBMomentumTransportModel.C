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

#include "incompressibleHFDIBMomentumTransportModel.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(incompressibleHFDIBMomentumTransportModel, 0);
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleHFDIBMomentumTransportModel::incompressibleHFDIBMomentumTransportModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    HFDIBmomentumTransportModel(U, alphaRhoPhi, phi, viscosity),
    alpha_(alpha),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::incompressibleHFDIBMomentumTransportModel>
Foam::incompressibleHFDIBMomentumTransportModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
{
    return HFDIBmomentumTransportModel::New<incompressibleHFDIBMomentumTransportModel>
    (
        geometricOneField(),
        geometricOneField(),
        U,
        phi,
        phi,
        viscosity
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::incompressibleHFDIBMomentumTransportModel::devSigma() const
{
    return devTau();
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleHFDIBMomentumTransportModel::divDevSigma(volVectorField& U) const
{
    return divDevTau(U);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::incompressibleHFDIBMomentumTransportModel::devTau() const
{
    NotImplemented;
    return devSigma();
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleHFDIBMomentumTransportModel::divDevTau
(
    volVectorField& U
) const
{
    NotImplemented;
    return divDevSigma(U);
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleHFDIBMomentumTransportModel::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;
    return divDevSigma(U);
}


// ************************************************************************* //
