
volScalarField& T = TSolid[i];
volScalarField& kappa = kappaSolid[i];

uniformDimensionedScalarField& alpha = alphaSolid[i];
uniformDimensionedScalarField& cpRef = cpRefSolid[i];
uniformDimensionedScalarField& rhoRef = rhoRefSolid[i];

kappa = alpha*cpRef*rhoRef;
kappa.correctBoundaryConditions();

fv::options& fvOptions = solidHeatSources[i];
