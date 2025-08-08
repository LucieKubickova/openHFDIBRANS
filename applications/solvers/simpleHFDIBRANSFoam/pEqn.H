volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

MRF.makeRelative(phiHbyA);
adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);

if (simple.consistent())
{
    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

UEqn.clear();

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

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
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                                 + fvc::div(phiFbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);
    
    pEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux() + phiFbyA;
    }
}

if (cutPhi)
{
    HFDIBRANS.cutPhiInBoundaryCells(phi);
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
U = HbyA - rAtU()*fvc::grad(p) + rAU*f;
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
