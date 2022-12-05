// compute force term
f = HFDIBRANS.surface()*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
f.correctBoundaryConditions();

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

    if (nCorr == 10) // should be done more generaly
    {
        phi = phiHbyA - pEqn.flux() + phiFbyA;
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

// update velocity -- should be done every time?... well yes, Ui-U should go to zero, it should not matter
U += 1.0*HFDIBRANS.surface()*(Ui - U);
//~ f += 1.0*surface*(Ui - U)/runTime.deltaT();

// Sync proc-proc boundaries
f.correctBoundaryConditions();
