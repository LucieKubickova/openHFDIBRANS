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

#include "HFDIBmomentumTransportModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallFvPatch.H"
#include "nearWallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HFDIBmomentumTransportModel, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

//~ Foam::typeIOobject<Foam::IOdictionary>
//~ Foam::HFDIBmomentumTransportModel::readModelDict
//~ (
    //~ const objectRegistry& obr,
    //~ const word& group,
    //~ bool registerObject
//~ )
//~ {
    //~ typeIOobject<IOdictionary> HFDIBmomentumTransport
    //~ (
        //~ IOobject::groupName(typeName, group),
        //~ obr.time().constant(),
        //~ obr,
        //~ IOobject::MUST_READ_IF_MODIFIED,
        //~ IOobject::NO_WRITE,
        //~ registerObject
    //~ );

    //~ if (HFDIBmomentumTransport.headerOk())
    //~ {
        //~ return HFDIBmomentumTransport;
    //~ }
    //~ else
    //~ {
        //~ typeIOobject<IOdictionary> turbulenceProperties
        //~ (
            //~ IOobject::groupName("turbulenceProperties", group),
            //~ obr.time().constant(),
            //~ obr,
            //~ IOobject::MUST_READ_IF_MODIFIED,
            //~ IOobject::NO_WRITE,
            //~ registerObject
        //~ );

        //~ if (turbulenceProperties.headerOk())
        //~ {
            //~ return turbulenceProperties;
        //~ }
        //~ else
        //~ {
            //~ return HFDIBmomentumTransport;
        //~ }
    //~ }
//~ }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HFDIBmomentumTransportModel::HFDIBmomentumTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
: momentumTransportModel(U, alphaRhoPhi, phi, viscosity)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::HFDIBmomentumTransportModel::phi() const
{
    return phi_;
}


const Foam::volScalarField::Boundary& Foam::HFDIBmomentumTransportModel::y() const
{
    return nearWallDist::New(mesh_).y();
}


bool Foam::HFDIBmomentumTransportModel::read()
{
    return regIOobject::read();
}


void Foam::HFDIBmomentumTransportModel::validate()
{}


void Foam::HFDIBmomentumTransportModel::correct()
{}


void Foam::HFDIBmomentumTransportModel::correct(openHFDIBRANS& HFDIBRANS)
{}


// ************************************************************************* //
