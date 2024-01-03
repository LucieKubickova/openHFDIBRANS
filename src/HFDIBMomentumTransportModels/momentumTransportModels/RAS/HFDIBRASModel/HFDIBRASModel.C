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

//~ #include "HFDIBRASModel.H"
#include "NewtonianViscosityModel.H"

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
    ),

    viscosityModel_
    (
        coeffDict_.found("viscosityModel")
      ? laminarModels::generalisedNewtonianViscosityModel::New
        (
            coeffDict_,
            viscosity,
            U
        )
      : autoPtr<laminarModels::generalisedNewtonianViscosityModel>
        (
            new laminarModels::generalisedNewtonianViscosityModels::Newtonian
            (
                coeffDict_,
                viscosity,
                U
            )
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

    const word modelType =
        modelDict.subDict("HFDIBRAS").lookupBackwardsCompatible<word>
        (
            {"model", "HFDIBRASModel"}
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
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, viscosity)
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
    viscosityModel_->correct();
    BasicHFDIBMomentumTransportModel::correct();
}


template<class BasicHFDIBMomentumTransportModel>
void Foam::HFDIBRASModel<BasicHFDIBMomentumTransportModel>::correct(openHFDIBRANS& HFDIBRANS)
{
    viscosityModel_->correct();
    BasicHFDIBMomentumTransportModel::correct(HFDIBRANS);
}

// ************************************************************************* //
