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
#include "wallDist.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBRASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBKOmegaSST<BasicMomentumTransportModel>::HFDIBKOmegaSST::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBKOmegaSST<BasicMomentumTransportModel>::HFDIBKOmegaSST::
F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBKOmegaSST<BasicMomentumTransportModel>::HFDIBKOmegaSST::
F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
HFDIBKOmegaSST<BasicMomentumTransportModel>::HFDIBKOmegaSST::
F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicMomentumTransportModel>
void HFDIBKOmegaSST<BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void HFDIBKOmegaSST<BasicMomentumTransportModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
HFDIBKOmegaSST<BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
HFDIBKOmegaSST<BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
HFDIBKOmegaSST<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
HFDIBKOmegaSST<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
HFDIBKOmegaSST<BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

template<class BasicMomentumTransportModel>
void HFDIBKOmegaSST<BasicMomentumTransportModel>::matrixManipulate
(
    fvScalarMatrix& eqn,
    volScalarField& phi,
    volScalarField& surface
)
{
    DynamicList<label> cells;
    DynamicList<scalar> phis;

    forAll(lambda_, cellI)
    {
        if (surface[cellI] == 1.0)
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
HFDIBKOmegaSST<BasicMomentumTransportModel>::HFDIBKOmegaSST
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

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
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
    kSurface_
    (
        IOobject
        (
            "HFDIBKOmegaSST::kSurface",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero",dimless,0.0)
    ),
    omegaSurface_
    (
        IOobject
        (
            "HFDIBKOmegaSST::omegaSurface",
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
    HFDIBRASDict.lookup("disSurfaceType") >> omegaSurfaceType_;
    HFDIBRASDict.lookup("kSurfaceType") >> kSurfaceType_;
    kBoundaryValue_ = readScalar(HFDIBRASDict.lookup("kBoundaryValue"));
    omegaBoundaryValue_ = readScalar(HFDIBRASDict.lookup("disBoundaryValue"));
    tolKEqn_ = readScalar(HFDIBRASDict.lookup("tolKEqn"));
    maxKEqnIters_ = readLabel(HFDIBRASDict.lookup("maxKEqnIters"));

    // bound
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool HFDIBKOmegaSST<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<HFDIBRASModel<BasicMomentumTransportModel>>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void HFDIBKOmegaSST<BasicMomentumTransportModel>::correct(openHFDIBRANS& HFDIBRANS)
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
    HFDIBRANS.createBaseSurface(kSurface_, kSurfaceType_, kBoundaryValue_);
    HFDIBRANS.createBaseSurface(omegaSurface_, omegaSurfaceType_, omegaBoundaryValue_);

    eddyViscosity<HFDIBRASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    volScalarField::Internal gamma(this->gamma(F1));
    volScalarField::Internal beta(this->beta(F1));

    // HFDIB: update uTau
    HFDIBRANS.updateUTau(k_);

    // HFDIB: correct omega and G
    HFDIBRANS.correctOmegaG(omega_, G, U, k_, nu_, omegaSurface_);

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
     ==
        alpha()*rho()*gamma
       *min
        (
            GbyNu,
            (c1_/a1_)*betaStar_*omega_()
           *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
        )
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
      - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
      - fvm::SuSp
        (
            alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
            omega_
        )
      + Qsas(S2(), gamma, beta)
      + omegaSource()
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

    // HFDIBRANS: matrix manipulate
    matrixManipulate(omegaEqn.ref(), omega_, omegaSurface_);

    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);

    // HFDIBRANS: compute imposed field for the turbulent kinetic energy
    HFDIBRANS.computeKi(k_, ki_, nu_);
            
    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, F23), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.relax();
    fvOptions.constrain(kEqn);

    for (label nCorr = 0; nCorr < maxKEqnIters_; nCorr++)
    {
        kQ_ = kSurface_*(kEqn.A()*ki_ - kEqn.H());
        solve(kEqn == kQ_);

        Info << "HFDIBRANS: Max error in k -> ki is " << (max(kSurface_*(ki - k)).value()) << endl;

        if (max(kSurface_*(ki_ - k_)).value() < tolKEqn_)
        {
            Info << "HFDIBRANS: k converged to ki within max tolerance " << tolKEqn_ << endl;
            break;
        }

        // apply correction
        k_ += 1.0*kSurface_*(ki_ - k_);
    }

    fvOptions.correct(k_);
    bound(k_, this->kMin_);
    
    correctNut(S2, F23);
    HFDIBRANS.correctNut(k_, nu_);
}

template<class BasicMomentumTransportModel>
void HFDIBKOmegaSST<BasicMomentumTransportModel>::correct()
{
    Info << "HFDIBKOmegaSST::correct() not implemented" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBRASModels
} // End namespace Foam

// ************************************************************************* //
