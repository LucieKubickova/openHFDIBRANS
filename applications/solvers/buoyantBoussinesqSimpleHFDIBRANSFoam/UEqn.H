// Solve the momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
    // NOTE: g -- nope
);

UEqn.relax();

fvOptions.constrain(UEqn);

HFDIBRANS.computeUi(U, Ui, surfaceTypeU);
HFDIBRANS.updateSurface(surfaceU, surfaceTypeU);

if (simple.momentumPredictor())
{
    for (label nCorr = 0; nCorr < maxEqnItersU; nCorr++)
    {
        f = surfaceU*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
          + f
        );

        Info << "HFDIBRANS: Max error in U -> Ui is " << (max(mag(surfaceU*(Ui - U))).value()) << endl;

        if (max(mag(surfaceU*(Ui - U))).value() < tolEqnU)
        {
            Info << "HFDIBRANS: U converged to Ui within max tolerance " << tolEqnU << endl;
            break;
        }

        // apply correction 
        U += 1.0*surfaceU*(Ui - U);
    }

    fvOptions.correct(U);
}
