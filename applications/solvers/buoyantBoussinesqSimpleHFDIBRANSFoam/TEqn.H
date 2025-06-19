// HFDIBRANS -- NOTE: this before or inside non-orthogonal loop?
HFDIBRANS.computeTi(T, Ti, surfaceTypeT, TIn);
HFDIBRANS.updateSurface(surfaceT, surfaceTypeT);

alphat = turbulence->nut()/Prt;
alphat.correctBoundaryConditions();

// HFDIBRANS
const tmp<volScalarField> tnu = turbulence->nu();
const volScalarField& nu = tnu();
HFDIBRANS.correctAlphat(alphat, nu, surfaceT);

volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

fvScalarMatrix TEqn
(
    fvm::div(phi, T)
  - fvm::laplacian(alphaEff, T)
 ==
    radiation->ST(rhoCpRef, T)
  + fvOptions(T)
);

TEqn.relax();

fvOptions.constrain(TEqn);

for (label nCorr = 0; nCorr < maxEqnItersT; nCorr++)
{
    // HFDIBRANS: update source
    Tq = surfaceT*(TEqn.A()*Ti - TEqn.H());
    solve(TEqn == Tq);

    Info << "HFDIB: Max error in T -> Ti is " << (max(mag(surfaceT*(Ti - T))).value()) << endl;

    if (max(surfaceT*(Ti - T)).value() < tolEqnT)
    {
        Info << "HFDIB: T converged to Ti within max tolerance " << tolEqnT << endl;
        break;
    }

    // apply correction
    T += 1.0*surfaceT*(Ti - T);
}

radiation->correct();

fvOptions.correct(T);

rhok = 1.0 - beta*(T - TRef);
