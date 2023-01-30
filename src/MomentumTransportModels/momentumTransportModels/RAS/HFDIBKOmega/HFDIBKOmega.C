/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "fvOptions.H"
#include "bound.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBRASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void HFDIBKOmega<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions(); // here nutWallFunctions are used
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> HFDIBKOmega<BasicMomentumTransportModel>::kSource() const
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
tmp<fvScalarMatrix> HFDIBKOmega<BasicMomentumTransportModel>::omegaSource() const
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
void HFDIBKOmega<BasicMomentumTransportModel>::matrixManipulate
(
    fvScalarMatrix& eqn,
    volScalarField& phi,
    scalar threshold,
    openHFDIBRANS& HFDIBRANS
)
{
    DynamicList<label> cells;
    DynamicList<scalar> phis;

    forAll(lambda_, cellI)
    {
        if (HFDIBRANS.outSurface()[cellI] == 1.0)
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
HFDIBKOmega<BasicMomentumTransportModel>::HFDIBKOmega
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

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.072
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
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
    nu_(this->nu())
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool HFDIBKOmega<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<HFDIBRASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void HFDIBKOmega<BasicMomentumTransportModel>::correct(openHFDIBRANS& HFDIBRANS)
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

    eddyViscosity<HFDIBRASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // HFDIB: correct omega and G
    HFDIBRANS.correctOmegaG(omega_, G, U, k_, nu_);

    // write G field
    if (this->runTime_.writeTime())
    {
        G.write();
    }

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_) // tohle je nula
      + fvm::div(alphaRhoPhi, omega_) // az na upravovane bunky stejne
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_) // bez korekce nut dobre v ramci vypoctu wall functions
     ==
        omegaSource() // tohle je nula
      + gamma_*alpha()*rho()*G*omega_()/k_() // tohle je dobre
      - fvm::Sp(beta_*alpha()*rho()*omega_(), omega_) // tohle je dobre v ramci vypoctu wall functions
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_) // tohle je stejny, pro stejny divU, juch
      + fvOptions(alpha, rho, omega_) // tohle je nula
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

    matrixManipulate(omegaEqn.ref(), omega_, 1e-4, HFDIBRANS);

    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);

    // HFDIB: compute imposed field for the turbulent kinetic energy
    HFDIBRANS.computeKi(k_, ki_, nu_);
            
    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(alpha, rho, k_) // nula
      + fvm::div(alphaRhoPhi, k_) // dobre az na posledni radu bunek logicky
      - fvm::laplacian(alpha*rho*DkEff(), k_) // dobre az na predposledni a posledni radu bunek
     ==
      kSource() // proste nula
      + alpha()*rho()*G // tohle vypada dobre
      - fvm::Sp(Cmu_*alpha()*rho()*omega_(), k_) // tohle je spravne v ramci spravnosti omega
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_) // tohle je spravne se spravnym divU
      + fvOptions(alpha, rho, k_) // nula
    );

    kEqn.relax();
    fvOptions.constrain(kEqn);

    for (label nCorr = 0; nCorr < HFDIBRANS.maxKEqnIters(); nCorr++)
    {
        kQ_ = HFDIBRANS.outSurface()*(kEqn.A()*ki_ - kEqn.H());
        solve(kEqn == kQ_);

        if (max(HFDIBRANS.outSurface()*(ki_ - k_)).value() < HFDIBRANS.tolKEqn())
        {
            Info << "HFDIBRAS: k converged to ki within max tolerance " << HFDIBRANS.tolKEqn() << endl;
            break;
        }

        // apply correction
        k_ += 1.0*HFDIBRANS.outSurface()*(ki_ - k_);
    }

    fvOptions.correct(k_);
    bound(k_, this->kMin_);
    
    correctNut(); // HFDIB: correct nut inside
}


template<class BasicMomentumTransportModel>
void HFDIBKOmega<BasicMomentumTransportModel>::correct()
{
    Info << "HFDIBKOmega::correct() not implemented" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBRASModels
} // End namespace Foam

// ************************************************************************* //

