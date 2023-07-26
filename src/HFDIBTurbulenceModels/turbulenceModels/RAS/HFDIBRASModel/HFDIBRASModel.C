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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::HFDIBRASModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::HFDIBRASModel<BasicTurbulenceModel>::HFDIBRASModel
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
    BasicTurbulenceModel
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
    HFDIBRASDict_(this->subOrEmptyDict("HFDIBRAS")),
    turbulence_(HFDIBRASDict_.getOrDefault<Switch>("turbulence", true)),
    printCoeffs_(HFDIBRASDict_.getOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(HFDIBRASDict_.optionalSubDict(type + "Coeffs")),
    kMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kMin",
            HFDIBRASDict_,
            sqr(dimVelocity),
            SMALL
        )
    ),
    epsilonMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "epsilonMin",
            HFDIBRASDict_,
            kMin_.dimensions()/dimTime,
            SMALL
        )
    ),
    omegaMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaMin",
            HFDIBRASDict_,
            dimless/dimTime,
            SMALL
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::HFDIBRASModel<BasicTurbulenceModel>>
Foam::HFDIBRASModel<BasicTurbulenceModel>::New
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

    const dictionary& dict = modelDict.subDict("HFDIBRAS");

    const word modelType
    (
        // "RASModel" for v2006 and earlier
        dict.getCompat<word>("model", {{"HFDIBRASModel", -2006}})
    );

    Info<< "Selecting HFDIBRAS turbulence model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "HFDIBRAS model",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<HFDIBRASModel>
    (
        ctorPtr(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::HFDIBRASModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        HFDIBRASDict_ <<= this->subDict("HFDIBRAS");
        HFDIBRASDict_.readEntry("turbulence", turbulence_);

        coeffDict_ <<= HFDIBRASDict_.optionalSubDict(type() + "Coeffs");

        kMin_.readIfPresent(HFDIBRASDict_);
        epsilonMin_.readIfPresent(HFDIBRASDict_);
        omegaMin_.readIfPresent(HFDIBRASDict_);

        return true;
    }

    return false;
}

template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBRASModel<BasicTurbulenceModel>::epsilon() const
{
    const scalar Cmu = 0.09;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        Cmu*this->k()*this->omega()
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::HFDIBRASModel<BasicTurbulenceModel>::omega() const
{
    const scalar betaStar = 0.09;
    const dimensionedScalar k0(sqr(dimLength/dimTime), SMALL);

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->mesh_.time().timeName(),
            this->mesh_
        ),
        this->epsilon()/(betaStar*(this->k() + k0))
    );
}

template<class BasicTurbulenceModel>
void Foam::HFDIBRASModel<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


template<class BasicTurbulenceModel>
void Foam::HFDIBRASModel<BasicTurbulenceModel>::correct(openHFDIBRANS& HFDIBRANS)
{
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
