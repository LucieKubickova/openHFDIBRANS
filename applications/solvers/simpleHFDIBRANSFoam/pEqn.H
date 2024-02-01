volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

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
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
U = HbyA - rAU*fvc::grad(p) + rAU*f;
U.correctBoundaryConditions();
fvOptions.correct(U);

// HFDIBRANS
if (cutVelocity)
{
    HFDIBRANS.cutUInBoundaryCells(U);
}

if (enforceVelocity)
{
    HFDIBRANS.enforceUiInBody(U, Ui);
}

// Sync proc-proc boundaries
f.correctBoundaryConditions();

// update new field
//~ f += 0.8*lambda*(Ui - U);
//~ f += surface*(Ui - U)/dummyTime;
