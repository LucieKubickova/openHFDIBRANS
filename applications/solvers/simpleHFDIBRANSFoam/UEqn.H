// Momentum predictor
fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    fvOptions(U)
  + g
);

UEqn.relax();

fvOptions.constrain(UEqn);

HFDIBRANS.computeUi(U, Ui, surfaceType);
Ui.correctBoundaryConditions(); // HERE NEW
HFDIBRANS.updateSurface(surface, surfaceType);
surface.correctBoundaryConditions(); // HERE NEW

if (simple.momentumPredictor())
{
    for (label nCorr = 0; nCorr < maxEqnIters; nCorr++)
    {
        f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        if (cutForce)
        {
            HFDIBRANS.cutFInBoundaryCells(f);
        }
        solve(UEqn == -fvc::grad(p) + f);

        Info << "HFDIBRANS: Max error in U -> Ui is " << (max(mag(surface*(Ui - U))).value()) << endl;

        if (max(mag(surface*(Ui - U))).value() < tolEqn)
        {
            Info << "HFDIBRANS: U converged to Ui within max tolerance " << tolEqn << endl;
            break;
        }

        // apply correction 
        U += 1.0*surface*(Ui - U);
    }

    fvOptions.correct(U);
}
