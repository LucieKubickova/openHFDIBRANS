
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

// calculate thermal diffusivity
alphaEff = turbulence->nu()/Pr + alphat;
HFDIBRANS.correctAlphaInSolid(alphaEff);
alphaEff.correctBoundaryConditions();

// interpolate alpha and correct
surfaceScalarField alphaSurface = fvc::interpolate(alphaEff);
HFDIBRANS.correctAlphaSurface(alphaSurface, alphaEff, rhoCp);

// calculate thermal conductivity
kappa = alphaEff*rhoCp;
kappa.correctBoundaryConditions();
surfaceScalarField kappaSurface = fvc::interpolate(kappa);
HFDIBRANS.correctKappaSurface(kappaSurface, kappa, kappa);

// interpolate rhoCp to surface
surfaceScalarField rhoCpSurface = linearInterpolate(rhoCp);
rhoCpSurface.correctBoundaryConditions();

surfaceScalarField lambdaSurface = fvc::interpolate(lambda);
HFDIBRANS.correctRhoCpSurface(rhoCpSurface, rhoCp);

fvScalarMatrix TEqn
(
    (1 - lambda)*fvm::div(rhoCpSurface*phi, T)
  - fvm::laplacian(kappaSurface, T)
  ==
    rhoCp*fvOptions(T)
);

//~ fvScalarMatrix TEqn
//~ (
//~     (1 - lambda)*fvm::div(phi, T)
//~   - fvm::laplacian(alphaSurface, T)
//~  ==
//~     radiation->ST(rhoCpRef, T)
//~   + fvOptions(T)
//~ );

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
