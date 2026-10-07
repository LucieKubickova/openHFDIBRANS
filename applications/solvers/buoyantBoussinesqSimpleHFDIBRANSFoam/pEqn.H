volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

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

// calculate lambda phi
surfaceScalarField phiLambda
(
    "phiLambda",
    phi*fvc::interpolate(lambda)
);

while (simple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
                                   - fvc::div(phiLambda)
    );

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

    p_rghEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA - p_rghEqn.flux() - phiLambda;

        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

#include "continuityErrs.H"

p = p_rgh + rho*gh;

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

// Sync proc-proc boundaries
//~ f.correctBoundaryConditions();
