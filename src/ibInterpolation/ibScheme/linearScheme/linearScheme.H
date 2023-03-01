/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License
openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with (level-set) DEM (Isoz and Sourek, 2020)

SourceFiles
    openHFDIBDEM.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef linearScheme_H
#define linearScheme_H

#include "ibScheme.H"

using namespace Foam;

//---------------------------------------------------------------------------//
class linearScheme
:
    public ibScheme
{
    public:

    // Constructors ----------------------------------------------------
    linearScheme();

    // Destructors -----------------------------------------------------
    ~linearScheme();

    // Public methods --------------------------------------------------
    // -- wrapper on interpolate functions
    scalar interpolate
    (
        volScalarField& phi,
        interpolation<scalar>& interpPhi,
        const volScalarField& body,
        scalar& dirichletVal,
        scalar& scale,
        scalar& ds,
        interpolationInfo& intInfo,
        label& cellI
    )
    {
        return interpolateT(phi, interpPhi, body, dirichletVal, scale, ds, intInfo, cellI);
    }

    vector interpolate
    (
        volVectorField& phi,
        interpolation<vector>& interpPhi,
        const volScalarField& body,
        vector& dirichletVal,
        scalar& scale,
        scalar& ds,
        interpolationInfo& intInfo,
        label& cellI
    )
    {
        return interpolateT(phi, interpPhi, body, dirichletVal, scale, ds, intInfo, cellI);
    }

    tensor interpolate
    (
        volTensorField& phi,
        interpolation<tensor>& interpPhi,
        const volScalarField& body,
        tensor& dirichletVal,
        scalar& scale,
        scalar& ds,
        interpolationInfo& intInfo,
        label& cellI
    )
    {
        return interpolateT(phi, interpPhi, body, dirichletVal, scale, ds, intInfo, cellI);
    }
    
    private:

    // -- templated interpolate funtion
    template <typename Type, typename volTypeField>
    Type interpolateT(volTypeField& phi, interpolation<Type>& interpPhi, const volScalarField& body, Type& dirichletVal, scalar& scale, scalar& ds, interpolationInfo& intInfo, label& cellI);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "linearSchemeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
