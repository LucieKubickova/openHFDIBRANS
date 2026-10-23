// Initialise fluid field pointer lists
//~ PtrList<rhoThermo> thermoFluid(fluidRegions.size());
//~ PtrList<volScalarField> rhoFluid(fluidRegions.size());
//~ PtrList<volVectorField> UFluid(fluidRegions.size());
//~ PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
//~ PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
//~ PtrList<volScalarField> ghFluid(fluidRegions.size());
//~ PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
//~ PtrList<volScalarField> p_rghFluid(fluidRegions.size());
//~ PtrList<radiation::radiationModel> radiation(fluidRegions.size());

PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulence(fluidRegions.size());
PtrList<singlePhaseTransportModel> laminarTransportFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volScalarField> kappaFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> cpRefFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> rhoRefFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> TRefFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> betaFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> PrFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> PrtFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size(), -1);
List<scalar> pRefValueFluid(fluidRegions.size(), Zero);

PtrList<dimensionedScalar> rhoMax(fluidRegions.size());
PtrList<dimensionedScalar> rhoMin(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

PtrList<fvVectorMatrix> UEqFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to TFluid\n" << endl;
    TFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to laminarTransportFluid\n" << endl;
    laminarTransportFluid.set
    (
        i,
        new singlePhaseTransportModel
        (
            UFluid[i],
            phiFluid[i]
        )
    );

    Info<< "    Adding to TRefFluid\n" << endl;
    TRefFluid.set
    (
        i, 
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "TRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("TRef", dimTemperature, laminarTransportFluid[i])
        )
    );

    Info<< "    Adding to cpRefFluid\n" << endl;
    cpRefFluid.set
    (
        i, 
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "cpRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("cpRef", dimLength*dimLength/dimTime/dimTime/dimTemperature, laminarTransportFluid[i])
        )
    );

    Info<< "    Adding to rhoRefFluid\n" << endl;
    rhoRefFluid.set
    (
        i, 
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "rhoRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("rhoRef", dimMass/dimLength/dimLength/dimLength, laminarTransportFluid[i])
        )
    );

    Info<< "    Adding to betaFluid\n" << endl;
    betaFluid.set
    (
        i, 
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "beta",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("beta", dimless/dimTemperature, laminarTransportFluid[i])
        )
    );

    Info<< "    Adding to PrFluid\n" << endl;
    PrFluid.set
    (
        i, 
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "Pr",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("Pr", dimless, laminarTransportFluid[i])
        )
    );

    Info<< "    Adding to PrtFluid\n" << endl;
    PrtFluid.set
    (
        i, 
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "Prt",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("Prt", dimless, laminarTransportFluid[i])
        )
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            1.0 - betaFluid[i]*(TFluid[i] - TRefFluid[i])
        )
    );

    Info<< "    Adding to kappaFluid\n" << endl;
    kappaFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "kappa",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            //~ dimensionedScalar("kappa", dimLength*dimLength*dimLength*dimLength/dimTime/dimTime/dimTime/dimTemperature, 1.0)
            dimensionedScalar("kappa", dimMass*dimLength/dimTime/dimTime/dimTime/dimTemperature, 1.0)
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i, 
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, laminarTransportFluid[i])
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );
    
    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            UFluid[i],
            phiFluid[i],
            laminarTransportFluid[i]
        ).ptr()
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            p_rghFluid[i] + rhoFluid[i]*ghFluid[i]
        )
    );

    // Force p_rgh to be consistent with p
    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    radiation.set
    (
        i,
        radiation::radiationModel::New(TFluid[i])
    );

    const dictionary& simpleDict =
        fluidRegions[i].solutionDict().subDict("SIMPLE");

    setRefCell
    (
        pFluid[i],
        p_rghFluid[i],
        simpleDict,
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    UEqFluid.set
    (
        i,
        new fvVectorMatrix(UFluid[i], dimForce)
    );

    turbulence[i].validate();
}
