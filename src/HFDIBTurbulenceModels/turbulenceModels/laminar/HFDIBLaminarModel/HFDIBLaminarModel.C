/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

template<class BasicHFDIBTurbulenceModel>
void Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicHFDIBTurbulenceModel>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::HFDIBLaminarModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicHFDIBTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    HFDIBLaminarDict_(this->subOrEmptyDict("laminar")), // Note: do not change to HFDIBLaminar
    printCoeffs_(HFDIBLaminarDict_.getOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(HFDIBLaminarDict_.optionalSubDict(type + "Coeffs"))
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicHFDIBTurbulenceModel>
Foam::autoPtr<Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    const IOdictionary modelDict
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false // Do not register
        )
    );

    const dictionary* dictptr = modelDict.findDict("laminar"); // Note: do not change to HFDIBLaminar

    if (dictptr)
    {
        const dictionary& dict = *dictptr;

        const word modelType
        (
            // HFDIBLaminarModel -> model (after v2006)
            dict.getCompat<word>("model", {{"HFDIBLaminarModel", -2006}})
        );

        Info<< "Selecting HFDIBLaminar stress model " << modelType << endl;

        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "HFDIBLaminar model",
                modelType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<HFDIBLaminarModel>
        (
            ctorPtr
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport, propertiesName)
        );
    }
    else
    {
        Info<< "Selecting HFDIBLaminar stress model "
            << HFDIBLaminarModels::HFDIBStokes<BasicHFDIBTurbulenceModel>::typeName << endl;

        return autoPtr<HFDIBLaminarModel>
        (
            new HFDIBLaminarModels::HFDIBStokes<BasicHFDIBTurbulenceModel>
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport,
                propertiesName
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicHFDIBTurbulenceModel>
bool Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::read()
{
    if (BasicHFDIBTurbulenceModel::read())
    {
        HFDIBLaminarDict_ <<= this->subDict("laminar");

        coeffDict_ <<= HFDIBLaminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }

    return false;
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::nut() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("nut", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimViscosity, Zero)
    );
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), Zero)
    );
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::nuEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("nuEff", this->alphaRhoPhi_.group()), this->nu()
        )
    );
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::k() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions())/dimTime, Zero)
    );
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::omega() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::R() const
{
    return tmp<volSymmTensorField>::New
    (
        IOobject
        (
            IOobject::groupName("R", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedSymmTensor(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicHFDIBTurbulenceModel>
void Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::correct()
{
    BasicHFDIBTurbulenceModel::correct();
}


template<class BasicHFDIBTurbulenceModel>
void Foam::HFDIBLaminarModel<BasicHFDIBTurbulenceModel>::correct(openHFDIBRANS& HFDIBRANS)
{
    BasicHFDIBTurbulenceModel::correct();
}

// ************************************************************************* //
