/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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

//~ #include "HFDIBLaminarModel.H"
#include "HFDIBStokes.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
void Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::printCoeffs
(
    const word& type
)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::HFDIBLaminarModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    BasicHFDIBMomentumTransportModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    HFDIBLaminarDict_(this->subOrEmptyDict("laminar")), // Note: do not change to HFDIBLaminar
    printCoeffs_(HFDIBLaminarDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(HFDIBLaminarDict_.optionalSubDict(type + "Coeffs"))
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
Foam::autoPtr<Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
{
    const IOdictionary modelDict
    (
        HFDIBmomentumTransportModel::readModelDict
        (
            U.db(),
            alphaRhoPhi.group()
        )
    );

    if (modelDict.found("laminar")) // Note: do not change to HFDIBLaminar
    {
        const word modelType =
            modelDict.subDict("laminar").lookupBackwardsCompatible<word> // Note: do not change to HFDIBLaminar
            (
                {"model", "laminarModel"} // Note: do not change to HFDIBLaminar
            );

        Info<< "Selecting HFDIBLaminar stress model " << modelType << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown HFDIBLaminarModel type "
                << modelType << nl << nl
                << "Valid HFDIBLaminarModel types:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<HFDIBLaminarModel>
        (
            cstrIter()
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                viscosity
            )
        );
    }
    else
    {
        Info<< "Selecting HFDIBLaminar stress model "
            << HFDIBLaminarModels::HFDIBStokes<BasicHFDIBMomentumTransportModel>::typeName
            << endl;

        return autoPtr<HFDIBLaminarModel>
        (
            new HFDIBLaminarModels::HFDIBStokes<BasicHFDIBMomentumTransportModel>
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                viscosity
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
bool Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::read()
{
    if (BasicHFDIBMomentumTransportModel::read())
    {
        HFDIBLaminarDict_ <<= this->subDict("laminar"); // Note: do not change to HFDIBLaminar

        coeffDict_ <<= HFDIBLaminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicHFDIBMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::nut() const
{
    return volScalarField::New
    (
        IOobject::groupName("nut", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimViscosity, 0)
    );
}


template<class BasicHFDIBMomentumTransportModel>
Foam::tmp<Foam::scalarField>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), 0.0)
    );
}


template<class BasicHFDIBMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::k() const
{
    return volScalarField::New
    (
        IOobject::groupName("k", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(sqr(dimVelocity), 0)
    );
}


template<class BasicHFDIBMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(sqr(dimVelocity)/dimTime, 0)
    );
}


template<class BasicHFDIBMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::omega() const
{
    return volScalarField::New
    (
        IOobject::groupName("omega", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimless/dimTime, 0)
    );
}


template<class BasicHFDIBMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::sigma() const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("sigma", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedSymmTensor(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicHFDIBMomentumTransportModel>
void Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::correct()
{
    BasicHFDIBMomentumTransportModel::correct();
}


template<class BasicHFDIBMomentumTransportModel>
void Foam::HFDIBLaminarModel<BasicHFDIBMomentumTransportModel>::correct(openHFDIBRANS& HFDIBRANS)
{
    BasicHFDIBMomentumTransportModel::correct();
}


// ************************************************************************* //
