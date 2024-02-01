// HFDIBRANS update
HFDIBRANS.computeUi(U, Ui);
HFDIBRANS.updateSurface(surface, surfaceType);

// Momentum predictor
for (label nCorr = 0; nCorr < maxUEqnIters; nCorr++)
{
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

    // compute LU
    scalarField LU = UEqn.lower() + UEqn.upper();
    volVectorField deltaU = Ui - U;
    volVectorField LUDeltaU
    (
        IOobject
        (
            "LUDeltaU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        deltaU/dummyTime
    );
    forAll(mesh.C(), cellI)
    {
        LUDeltaU[cellI] *= LU[cellI];
    }

    // compute force
    //~ f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
    f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p) - LUDeltaU);
    //~ f += surface*(Ui - U)/dummyTime;
    //~ f = surface*(Ui - U)/dummyTime;
    //~ HFDIBRANS.cutFInBoundaryCells(f);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p) + f);

        if (max(mag(surface*(Ui - U))).value() < tolUEqn)
        {
            Info << "HFDIBRAS: U converged to Ui within max tolerance " << tolUEqn << endl;
            break;
        }

        // apply correction 
        U += 1.0*surface*(Ui - U);
    }
}

fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    fvOptions(U)
  + g
);

fvOptions.correct(U);
