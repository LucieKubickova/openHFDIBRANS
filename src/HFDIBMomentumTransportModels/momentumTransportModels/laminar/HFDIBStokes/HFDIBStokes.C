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

template<class BasicMomentumTransportModel>
HFDIBStokes<BasicMomentumTransportModel>::HFDIBStokes
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
:
    linearViscousStress<HFDIBLaminarModel<BasicMomentumTransportModel>>
    (
        typeName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
const dictionary&
HFDIBStokes<BasicMomentumTransportModel>::coeffDict() const
{
    return dictionary::null;
}


template<class BasicMomentumTransportModel>
bool HFDIBStokes<BasicMomentumTransportModel>::read()
{
    return true;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBStokes<BasicMomentumTransportModel>::nut() const
{
    return volScalarField::New
    (
        IOobject::groupName("nut", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimViscosity, 0)
    );
}


template<class BasicMomentumTransportModel>
tmp<scalarField>
HFDIBStokes<BasicMomentumTransportModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), 0.0)
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBStokes<BasicMomentumTransportModel>::nuEff() const
{
    return volScalarField::New
    (
        IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
        this->nu()
    );
}


template<class BasicMomentumTransportModel>
tmp<scalarField>
HFDIBStokes<BasicMomentumTransportModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBStokes<BasicMomentumTransportModel>::k() const
{
    return volScalarField::New
    (
        IOobject::groupName("k", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions()), 0)
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBStokes<BasicMomentumTransportModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions())/dimTime, 0)
    );
}


template<class BasicMomentumTransportModel>
tmp<volSymmTensorField>
HFDIBStokes<BasicMomentumTransportModel>::sigma() const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("R", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedSymmTensor(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicMomentumTransportModel>
void HFDIBStokes<BasicMomentumTransportModel>::correct()
{
    HFDIBLaminarModel<BasicMomentumTransportModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBLaminarModels
} // End namespace Foam

// ************************************************************************* //
