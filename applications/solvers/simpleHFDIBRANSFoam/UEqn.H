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
Ui.correctBoundaryConditions();
HFDIBRANS.updateSurface(surface, surfaceType);
surface.correctBoundaryConditions();
HFDIBRANS.correctSurfaceByNormal(normSurface, surface, normCorrLimit);
normSurface.correctBoundaryConditions();

if (simple.momentumPredictor())
{
    for (label nCorr = 0; nCorr < maxEqnIters; nCorr++)
    {
        if (useNormSurface)
        {
            f = normSurface & (UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        }

        else
        {
            f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        }

        if (cutForce)
        {
            HFDIBRANS.cutFInBoundaryCells(f);
        }
        solve(UEqn == -fvc::grad(p) + f);

        scalar error(0.0);
        if (useNormSurface)
        {
            error = max(mag(normSurface & (Ui - U))).value();
        }

        else
        {
            error = max(mag(surface*(Ui - U))).value();
        }
        Info << "HFDIBRANS: Max error in U -> Ui is " << error << endl;

        if (error < tolEqn)
        {
            Info << "HFDIBRANS: U converged to Ui within max tolerance " << tolEqn << endl;
            break;
        }

        // apply correction 
        if (useNormSurface)
        {
            U += 1.0*normSurface & (Ui - U);
        }

        else
        {
            U += 1.0*surface*(Ui - U);
        }
    }

    fvOptions.correct(U);
}
