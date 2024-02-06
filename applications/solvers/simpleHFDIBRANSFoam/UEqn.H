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

HFDIBRANS.computeUi(U, Ui);
HFDIBRANS.updateSurface(surface, surfaceType);

if (simple.momentumPredictor())
{
    for (label nCorr = 0; nCorr < maxUEqnIters; nCorr++)
    {
        f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        solve(UEqn == -fvc::grad(p) + f);

        Info << "HFDIBRANS: Max error in U -> Ui is " << (max(mag(surface*(Ui - U))).value()) << endl;

        if (max(mag(surface*(Ui - U))).value() < tolUEqn)
        {
            Info << "HFDIBRANS: U converged to Ui within max tolerance " << tolUEqn << endl;
            break;
        }

        // apply correction 
        U += 1.0*surface*(Ui - U);
    }

    fvOptions.correct(U);
}
