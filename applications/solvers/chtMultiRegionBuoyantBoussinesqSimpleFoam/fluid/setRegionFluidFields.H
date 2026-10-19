
    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& T = TFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    incompressible::turbulenceModel& turb = turbulence[i];
    singlePhaseTransportModel& lamTrans = laminarTransportFluid[i];

    volScalarField& p = pFluid[i];
    volScalarField& p_rgh = p_rghFluid[i];

    volScalarField& kappa = kappaFluid[i];

    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    uniformDimensionedScalarField& beta = betaFluid[i];
    uniformDimensionedScalarField& TRef = TRefFluid[i];
    uniformDimensionedScalarField& cpRef = cpRefFluid[i];
    uniformDimensionedScalarField& rhoRef = rhoRefFluid[i];
    uniformDimensionedScalarField& Pr = PrFluid[i];
    uniformDimensionedScalarField& Prt = PrtFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];
