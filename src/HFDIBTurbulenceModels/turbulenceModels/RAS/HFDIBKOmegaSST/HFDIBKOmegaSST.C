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

template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::HFDIBKOmegaSST::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    // HFDIB: correct value
    tmp<volScalarField> y = max // NOTE: wrong
    (
        y_,
        yIB_
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y()),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y())*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y()))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::HFDIBKOmegaSST::F2() const
{
    // HFDIB: correct value
    tmp<volScalarField> y = max // NOTE: wrong
    (
        y_,
        yIB_
    );

    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y()),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y())*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::HFDIBKOmegaSST::F3() const
{
    // HFDIB: correct value
    tmp<volScalarField> y = max // NOTE: wrong
    (
        y_,
        yIB_
    );

    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y())),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::HFDIBKOmegaSST::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicHFDIBTurbulenceModel>
void HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::correctNut
(
    const volScalarField& S2
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicHFDIBTurbulenceModel>
void HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicHFDIBTurbulenceModel>
Foam::tmp<Foam::volScalarField> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::S2
(
    const volTensorField& gradU
) const
{
    return 2*magSqr(symm(gradU));
}


template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField::Internal> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField::Internal> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::epsilonByk
(
    const volScalarField& /* F1 not used */,
    const volTensorField& /* gradU not used */
) const
{
    return betaStar_*omega_();
}


template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField::Internal> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::GbyNu0
(
    const volTensorField& gradU,
    const volScalarField& /* S2 not used */
) const
{
    return tmp<volScalarField::Internal>::New
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        gradU() && dev(twoSymm(gradU()))
    );
}


template<class BasicHFDIBTurbulenceModel>
tmp<volScalarField::Internal> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega_()*max(a1_*omega_(), b1_*F2*sqrt(S2))
    );
}


template<class BasicHFDIBTurbulenceModel>
tmp<fvScalarMatrix> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        k_,
        dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
    );
}


template<class BasicHFDIBTurbulenceModel>
tmp<fvScalarMatrix> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        omega_,
        dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
    );
}


template<class BasicHFDIBTurbulenceModel>
tmp<fvScalarMatrix> HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>::New
    (
        omega_,
        dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
    );
}

template<class BasicTurbulenceModel>
void HFDIBKOmegaSST<BasicTurbulenceModel>::matrixManipulate
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

template<class BasicHFDIBTurbulenceModel>
HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::HFDIBKOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<HFDIBRASModel<BasicHFDIBTurbulenceModel>>
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

    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),
    yIB_
    (
        IOobject
        (
            "HFDIBKOmegaSST::yIB",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero",y_.dimensions(),0.0)
    ),
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

    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
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
    omegaGSurface_
    (
        IOobject
        (
            "HFDIBKOmegaSST::omegaGSurface",
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
    HFDIBRASDict.readEntry("kSurfaceType", kSurfaceType_);
    HFDIBRASDict.readEntry("disGSurfaceType", omegaGSurfaceType_);
    HFDIBRASDict.readEntry("kBoundaryValue", kBoundaryValue_);
    HFDIBRASDict.readEntry("disGBoundaryValue", omegaGBoundaryValue_);
    HFDIBRASDict.readEntry("tolKEqn", tolKEqn_);
    HFDIBRASDict.readEntry("maxKEqnIters", maxKEqnIters_);
    useKQ_ = HFDIBRASDict.lookupOrDefault<bool>("useKSource", true);
    correctFs_ = HFDIBRASDict.lookupOrDefault<bool>("correctFs", true);

    // bound
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicHFDIBTurbulenceModel>
void HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}


template<class BasicHFDIBTurbulenceModel>
bool HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::read()
{
    if (BasicHFDIBTurbulenceModel::read())
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

        setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicHFDIBTurbulenceModel>
void HFDIBKOmegaSST<BasicHFDIBTurbulenceModel>::correct(openHFDIBRANS& HFDIBRANS)
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

    BasicHFDIBTurbulenceModel::correct();

    // HFDIBRANS references
    HFDIBRANS.createBaseSurface(kSurface_, kSurfaceType_, kBoundaryValue_);
    HFDIBRANS.createBaseSurface(omegaGSurface_, omegaGSurfaceType_, omegaGBoundaryValue_);

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));
    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    // HFDIB: correct y_
    if (correctFs_)
    {
        HFDIBRANS.correctY(yIB_);
    }

    const volScalarField F1(this->F1(CDkOmega));
    const volScalarField F23(this->F23());
    //~ F1.rename("F1"); // Note (LK): not working for v2412
    //~ F23.rename("F23"); // Note (LK): not working for v2412

    if (this->runTime_.writeTime())
    {
        F1.write();
        F23.write();
        y_.write();
    }

    // HFDIB: correct F1 and F23
    //~ HFDIBRANS.correctF1(F1);
    //~ HFDIBRANS.correctF23(F23);

    {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = GbyNu(GbyNu0, F23(), S2());

        // HFDIB: update uTau
        HFDIBRANS.updateUTau(k_);

        // HFDIB: correct omega and G
        HFDIBRANS.correctOmegaG(omega_, G, U, k_, nu_, omegaGSurface_);

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

        // HFDIBRANS: matrix manipulate
        matrixManipulate(omegaEqn.ref(), omega_, omegaGSurface_);

        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
        //~ HFDIBRANS.bound(omega_, this->omegaMin_);
    }

    // HFDIBRANS: compute imposed field for the turbulent kinetic energy
    HFDIBRANS.computeKi(k_, ki_, nu_);
            
    {
        // Turbulent kinetic energy equation
        fvScalarMatrix kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + fvm::div(alphaRhoPhi, k_)
          - fvm::laplacian(alpha*rho*DkEff(F1), k_)
         ==
            alpha()*rho()*Pk(G)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
          - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
          + alpha()*rho()*betaStar_*omegaInf_*kInf_
          + kSource()
          + fvOptions(alpha, rho, k_)
        );

        tgradU.clear();

        kEqn.relax();
        fvOptions.constrain(kEqn);

        if (useKQ_)
        {
            for (label nCorr = 0; nCorr < maxKEqnIters_; nCorr++)
            {
                kQ_ = kSurface_*(kEqn.A()*ki_ - kEqn.H());
                solve(kEqn == kQ_);

                Info << "HFDIBRANS: Max error in k -> ki is " << (max(kSurface_*(ki_ - k_)).value()) << endl;

                if (max(kSurface_*(ki_ - k_)).value() < tolKEqn_)
                {
                    Info << "HFDIBRANS: k converged to ki within max tolerance " << tolKEqn_ << endl;
                    break;
                }

                // apply correction
                k_ += 1.0*kSurface_*(ki_ - k_);
            }

        }

        else
        {
            solve(kEqn);
        }
    }

    fvOptions.correct(k_);
    bound(k_, this->kMin_);
    //~ HFDIBRANS.bound(k_, this->kMin_);

    correctNut(S2);
    HFDIBRANS.correctNut(k_, nu_);
}

template<class BasicTurbulenceModel>
void HFDIBKOmegaSST<BasicTurbulenceModel>::correct()
{
    Info << "HFDIBKOmegaSST::correct() not implemented" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBRASModels
} // End namespace Foam

// ************************************************************************* //
