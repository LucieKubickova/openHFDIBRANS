/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "HFDIBStokes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBLaminarModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
HFDIBStokes<BasicHFDIBMomentumTransportModel>::HFDIBStokes
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    linearViscousStress<HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>>
    (
        typeName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
const dictionary& HFDIBStokes<BasicHFDIBMomentumTransportModel>::coeffDict() const
{
    return dictionary::null;
}


template<class BasicHFDIBMomentumTransportModel>
bool HFDIBStokes<BasicHFDIBMomentumTransportModel>::read()
{
    return true;
}


template<class BasicHFDIBMomentumTransportModel>
tmp<volScalarField> HFDIBStokes<BasicHFDIBMomentumTransportModel>::nuEff() const
{
    return volScalarField::New
    (
        IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
        this->nu()
    );
}


template<class BasicHFDIBMomentumTransportModel>
tmp<scalarField> HFDIBStokes<BasicHFDIBMomentumTransportModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicHFDIBMomentumTransportModel>
void HFDIBStokes<BasicHFDIBMomentumTransportModel>::correct()
{
    HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBLaminarModels
} // End namespace Foam

// ************************************************************************* //
