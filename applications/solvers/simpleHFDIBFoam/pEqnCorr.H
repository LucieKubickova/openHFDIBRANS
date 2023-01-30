for (label nCorr = 0; nCorr < HFDIBRANS.maxPEqnIters(); nCorr++)
{
    // auxiliar boolean
    bool toEnd(false);

    // compute force term
    //~ f = HFDIBRANS.outSurface()*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
    //~ f = HFDIBRANS.inSurface()*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
    //~ f.correctBoundaryConditions();
    
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);
    
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAU, MRF);
    
    volVectorField FbyA(rAU*f);
    surfaceScalarField phiFbyA
    (
        "phiFbyA",
        fvc::flux(FbyA)
    );
    
    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) ==  fvc::div(phiHbyA)
                                     + fvc::div(phiFbyA)
        );
        
        pEqn.setReference(pRefCell, pRefValue);
        
        pEqn.solve();

        //~ if (simple.finalNonOrthogonalIter())
        //~ {
            //~ phi = phiHbyA - pEqn.flux() + phiFbyA;

            // correct phi at the immersed boundary
            //~ HFDIBRANS.correctPhi(phi);
        //~ }
    //~ }

    //~ // check error in pressure boundary condition
    //~ if (HFDIBRANS.maxErrorInPBC(p) < 1e-4)
    //~ {
        //~ Info << "HFDIBRANS: p converged to zeroGradient boundary condition at the IB" << endl;
        //~ toEnd = true;
    //~ }

    //~ // correct pressure
    //~ HFDIBRANS.correctP(p);

    if (toEnd or nCorr == HFDIBRANS.maxPEqnIters()-1)
    {
        phi = phiHbyA - pEqn.flux() + phiFbyA;

        // correct phi at the immersed boundary
        //~ HFDIBRANS.correctPhi(phi);

        #include "continuityErrs.H"
        
        // Explicitly relax pressure for momentum corrector
        p.relax();
        // Explicitly relax the forcing term for momentum corrector
        //~ f.relax(); // what about this?
        
        // Momentum corrector
        U = HbyA - rAU*fvc::grad(p) + rAU*f;
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        
        // Sync proc-proc boundaries
        f.correctBoundaryConditions();
    }
    }
}
