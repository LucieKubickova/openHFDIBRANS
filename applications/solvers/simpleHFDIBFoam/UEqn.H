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
    for (label nCorr = 0; nCorr < HFDIBRANS.maxUEqnIters(); nCorr++)
    {
        f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        solve(UEqn == -fvc::grad(p) + f);

        if (max(mag(surface*(Ui - U))).value() < HFDIBRANS.tolUEqn())
        {
            Info << "HFDIBRAS: U converged to Ui within max tolerance " << HFDIBRANS.tolUEqn() << endl;
            break;
        }

        // apply correction 
        U += 1.0*surface*(Ui - U);
    }

    fvOptions.correct(U);
}
