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

#include "fvOptions.H"
#include "bound.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBRASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> HFDIBRealizableKE<BasicMomentumTransportModel>::rCmu
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    tmp<volSymmTensorField> tS = dev(symm(gradU));
    const volSymmTensorField& S = tS();

    volScalarField W
    (
        (2*sqrt(2.0))*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar(dimensionSet(0, 0, -3, 0, 0), small)
        )
    );

    tS.clear();

    volScalarField phis
    (
        (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)))
    );
    volScalarField As(sqrt(6.0)*cos(phis));
    volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

    return 1.0/(A0_ + As*Us*k_/epsilon_);
}


template<class BasicMomentumTransportModel>
void HFDIBRealizableKE<BasicMomentumTransportModel>::correctNut
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    this->nut_ = rCmu(gradU, S2, magS)*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicMomentumTransportModel>
void HFDIBRealizableKE<BasicMomentumTransportModel>::correctNut()
{
    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    volScalarField S2(modelName("S2"), 2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(modelName("magS"), sqrt(S2));

    correctNut(tgradU(), S2, magS);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> HFDIBRealizableKE<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
HFDIBRealizableKE<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}

template<class BasicMomentumTransportModel>
void HFDIBRealizableKE<BasicMomentumTransportModel>::matrixManipulate
(
    fvScalarMatrix& eqn,
    volScalarField& phi,
    scalar threshold
)
{
    DynamicList<label> cells;
    DynamicList<scalar> phis;

    forAll(lambda_, cellI)
    {
        if (surface_[cellI] == 1.0)
        {
            cells.append(cellI);
            phis.append(phi[cellI]);
        }
    }

    labelUList UCells(cells);
    UList<scalar> UPhis(phis);

    eqn.setValues(UCells, UPhis);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
HFDIBRealizableKE<BasicMomentumTransportModel>::HFDIBRealizableKE
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<HFDIBRASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.0
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.9
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.2
        )
    ),
    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_
    ),
    surface_
    (
        IOobject
        (
            "HFDIBRealizableKE::surface",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero",dimless,0.0)
    ),
    ki_
    (
        IOobject
        (
            "ki",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero",k_.dimensions(),0)
    ),
    kQ_
    (
        IOobject
        (
            "kQ",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    nu_(this->nu()),
    HFDIBDEMDict_
    (
        IOobject
        (
            "HFDIBDEMDict",
            this->runTime_.constant(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    // read dictionaries
    dictionary HFDIBRASDict = this->HFDIBRASDict_;
    HFDIBRASDict.lookup("surfaceType") >> surfaceType_;
    boundaryValue_ = readScalar(HFDIBRASDict.lookup("boundaryValue"));
    tolKEqn_ = readScalar(HFDIBRASDict.lookup("tolKEqn"));
    maxKEqnIters_ = readLabel(HFDIBRASDict.lookup("maxKEqnIters"));

    // bound
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool HFDIBRealizableKE<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<HFDIBRASModel<BasicMomentumTransportModel>>::read())
    {
        A0_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void HFDIBRealizableKE<BasicMomentumTransportModel>::correct(openHFDIBRANS& HFDIBRANS)
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // HFDIBRANS references
    HFDIBRANS.createBaseSurface(surface_, surfaceType_, boundaryValue_);

    eddyViscosity<HFDIBRASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        modelName("divU"),
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(modelName("S2"), 2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(modelName("magS"), sqrt(S2));

    volScalarField::Internal eta(modelName("eta"), magS()*k_()/epsilon_());
    volScalarField::Internal C1
    (
        modelName("C1"),
        max(eta/(scalar(5) + eta), scalar(0.43))
    );

    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    
    // HFDIB: correct epsilon and G
    HFDIBRANS.correctEpsilonG(epsilon_, G, U, k_, nu_, surface_);

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1*alpha()*rho()*magS()*epsilon_()
      - fvm::Sp
        (
            C2_*alpha()*rho()*epsilon_()/(k_() + sqrt(this->nu()()*epsilon_())),
            epsilon_
        )
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());

    // HFDIBRANS: matrix manipulate
    matrixManipulate(epsEqn.ref(), epsilon_, 1e-4);

    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // HFDIBRANS: compute imposed field for the turbulent kinetic energy
    HFDIBRANS.computeKi(k_, ki_, nu_);

    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.relax();
    fvOptions.constrain(kEqn);

    for (label nCorr = 0; nCorr < maxKEqnIters_; nCorr++)
    {
        kQ_ = surface_*(kEqn.A()*ki_ - kEqn.H());
        solve(kEqn == kQ_);

        if (max(surface_*(ki_ - k_)).value() < tolKEqn_)
        {
            Info << "HFDIBRAS: k converged to ki within max tolerance " << tolKEqn_ << endl;
            break;
        }

        // apply correction
        k_ += 1.0*surface_*(ki_ - k_);
    }
    
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(tgradU(), S2, magS);
    HFDIBRANS.correctNut(k_, nu_);
}

template<class BasicMomentumTransportModel>
void HFDIBRealizableKE<BasicMomentumTransportModel>::correct()
{
    Info << "HFDIBRealizableKE::correct() not implemented" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBRASModels
} // End namespace Foam

// ************************************************************************* //
