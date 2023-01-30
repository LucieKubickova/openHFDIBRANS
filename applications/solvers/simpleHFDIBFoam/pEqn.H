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

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux() + phiFbyA;

        // correct phi at the immersed boundary
        //~ HFDIBRANS.correctPhi(phi);
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();
// Explicitly relax the forcing term for momentum corrector
//~ f.relax(); // what about this?

// Momentum corrector
U = HbyA - rAU*fvc::grad(p) + rAU*f;
U.correctBoundaryConditions();
fvOptions.correct(U);

//~ U += 1.0*surface*(Ui - U);
//~ f += 1.0*surface*(Ui - U)/runTime.deltaT();
// Sync proc-proc boundaries
//~ f.correctBoundaryConditions();
