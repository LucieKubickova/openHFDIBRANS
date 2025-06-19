volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

MRF.makeRelative(phiHbyA);
adjustPhi(phiHbyA, U, p_rgh);

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

volVectorField FbyA(rAU*f);
surfaceScalarField phiFbyA
(
    "phiFbyA",
    fvc::flux(FbyA)
);

while (simple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
                                   + fvc::div(phiFbyA)
    );

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

    p_rghEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA - p_rghEqn.flux() + phiFbyA;

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf) + rAU*f;
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

#include "continuityErrs.H"

p = p_rgh + rhok*gh;

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rhok*gh;
}

// Sync proc-proc boundaries
f.correctBoundaryConditions();
