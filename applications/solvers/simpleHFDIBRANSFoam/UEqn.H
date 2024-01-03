    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevSigma(U)
     ==
        fvModels.source(U)
      + g
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    HFDIBRANS.computeUi(U, Ui);
    HFDIBRANS.updateSurface(surface, surfaceType);
    
    if (simple.momentumPredictor())
    {
        for (label nCorr = 0; nCorr < maxUEqnIters; nCorr++)
        {
            f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
            solve(UEqn == -fvc::grad(p) + f);
    
            if (max(mag(surface*(Ui - U))).value() < tolUEqn)
            {
                Info << "HFDIBRAS: U converged to Ui within max tolerance " << tolUEqn << endl;
                break;
            }
    
            // apply correction 
            U += 1.0*surface*(Ui - U);
        }
    
        fvConstraints.constrain(U);
    }
