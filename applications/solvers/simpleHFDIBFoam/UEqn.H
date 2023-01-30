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

HFDIBRANS.computeUi(U, Ui, turbulence->k(), turbulence->nu()); // problem for laminar, or... is it?

if (simple.momentumPredictor())
{
    for (label nCorr = 0; nCorr < HFDIBRANS.maxUEqnIters(); nCorr++)
    {
        //~ f = HFDIBRANS.outSurface()*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        f = HFDIBRANS.inSurface()*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        solve(UEqn == -fvc::grad(p) + f);

        //~ if (max(mag(HFDIBRANS.outSurface()*(Ui - U))).value() < HFDIBRANS.tolUEqn())
        if (max(mag(HFDIBRANS.inSurface()*(Ui - U))).value() < HFDIBRANS.tolUEqn())
        {
            Info << "HFDIBRAS: U converged to Ui within max tolerance " << HFDIBRANS.tolUEqn() << endl;
            break;
        }

        // apply correction 
        //~ U += 1.0*HFDIBRANS.outSurface()*(Ui - U);
        U += 1.0*HFDIBRANS.inSurface()*(Ui - U);
    }

    fvOptions.correct(U);
}
