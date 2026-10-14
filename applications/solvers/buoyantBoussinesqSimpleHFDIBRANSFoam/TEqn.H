
// HFDIBRANS: set up imposed temperature
HFDIBRANS.computeTi(T, Ti, surfaceTypeT, TIn, boundaryValT);
Ti.correctBoundaryConditions();

// HFDIBRANS: set up surface
HFDIBRANS.updateSurface(surfaceT, surfaceTypeT, boundaryValT);
surfaceT.correctBoundaryConditions();

alphat = turbulence->nut()/Prt;
alphat.correctBoundaryConditions();

// HFDIBRANS: correct alphat
const tmp<volScalarField> tnu = turbulence->nu();
const volScalarField& nu = tnu();
HFDIBRANS.correctAlphat(alphat, nu, surfaceT);

//~ volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);
alphaEff = turbulence->nu()/Pr + alphat;
HFDIBRANS.correctAlphaInSolid(alphaEff);
alphaEff.correctBoundaryConditions();

fvScalarMatrix TEqn
(
    (1 - lambda)*fvm::div(phi, T)
  - fvm::laplacian(alphaEff, T)
 ==
    radiation->ST(rhoCpRef, T)
  + fvOptions(T)
);

TEqn.relax();

fvOptions.constrain(TEqn);

// HFDIBRANS: immersed boundary induced source term
if (!solveTInSolid)
{
    Tq = surfaceT*(TEqn.A()*Ti - TEqn.H());
}

// HFDIBRANS: add to T equation
TEqn -= Tq;

solve(TEqn == Tq);

// calculate error
scalar errorT(0.0);
errorT = max(mag(surfaceT*(Ti - T))).value();
Info << "HFDIBRANS: Max error in T -> Ti is " << errorT << endl;

// apply correction
T += 1.0*surfaceT*(Ti - T);

radiation->correct();

fvOptions.correct(T);

rho = 1.0 - beta*(T - TRef);
