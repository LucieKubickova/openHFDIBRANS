// Solve the momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    fvOptions(U)
    // NOTE: g -- nope
);

UEqn.relax();

fvOptions.constrain(UEqn);

// HFDIBRANS: set up imposed velocity
HFDIBRANS.computeUi(U, Ui, surfaceTypeU, boundaryValU);
Ui.correctBoundaryConditions();

// HFDIBRANS: set up surfaces
HFDIBRANS.updateSurface(surfaceU, surfaceTypeU, boundaryValU);
surfaceU.correctBoundaryConditions();

// HFDIBRANS: immersed boundary induced source term
f = surfaceU*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));

// HFDIBRANS: add to U equation
UEqn -= f;

if (simple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    // Note (LK): following reconstruct needs fixing of ghf
    // and check of snGrad and magSf at immersed boundary
    //~ solve
    //~ (
    //~     UEqn
    //~   ==
    //~     fvc::reconstruct
    //~     (
    //~         (
    //~           - ghf*fvc::snGrad(rho)
    //~           - fvc::snGrad(p_rgh)
    //~         )*mesh.magSf()
    //~     )
    //~ );

    fvOptions.correct(U);
}

// calculate error
scalar error(0.0);
error = max(mag(surfaceU*(Ui - U))).value();
Info << "HFDIBRANS: Max error in U -> Ui is " << error << endl;

// apply correction
U += 1.0*surfaceU*(Ui - U);
