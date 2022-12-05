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

//~ #include "HFDIBRASModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
void Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>::HFDIBRASModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
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
        transport
    ),

    HFDIBRASDict_(this->subOrEmptyDict("HFDIBRAS")),
    turbulence_(HFDIBRASDict_.lookup("turbulence")),
    printCoeffs_(HFDIBRASDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(HFDIBRASDict_.optionalSubDict(type + "Coeffs")),

    kMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kMin",
            HFDIBRASDict_,
            sqr(dimVelocity),
            small
        )
    ),

    epsilonMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "epsilonMin",
            HFDIBRASDict_,
            kMin_.dimensions()/dimTime,
            small
        )
    ),

    omegaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaMin",
            HFDIBRASDict_,
            dimless/dimTime,
            small
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
Foam::autoPtr<Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>>
Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
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

    const word modelType
    (
        modelDict.subDict("HFDIBRAS").found("model")
      ? modelDict.subDict("HFDIBRAS").lookup("model")
      : modelDict.subDict("HFDIBRAS").lookup("HFDIBRASModel")
    );

    Info<< "Selecting HFDIBRAS turbulence model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown HFDIBRASModel type "
            << modelType << nl << nl
            << "Valid HFDIBRASModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<HFDIBRASModel>
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, transport)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicHFDIBMomentumTransportModel>
bool Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>::read()
{
    if (BasicHFDIBMomentumTransportModel::read())
    {
        HFDIBRASDict_ <<= this->subDict("HFDIBRAS");
        HFDIBRASDict_.lookup("turbulence") >> turbulence_;

        coeffDict_ <<= HFDIBRASDict_.optionalSubDict(type() + "Coeffs");

        kMin_.readIfPresent(HFDIBRASDict_);
        epsilonMin_.readIfPresent(HFDIBRASDict_);
        omegaMin_.readIfPresent(HFDIBRASDict_);

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicHFDIBMomentumTransportModel>
void Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>::correct()
{
    BasicHFDIBMomentumTransportModel::correct();
}


template<class BasicHFDIBMomentumTransportModel>
void Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>::correct(openHFDIBRANS& HFDIBRANS)
{
    correct();
}


// ************************************************************************* //
