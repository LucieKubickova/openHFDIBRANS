// Initialise fluid field pointer lists
//~ PtrList<rhoThermo> thermoFluid(fluidRegions.size());
//~ PtrList<volScalarField> rhoFluid(fluidRegions.size());
//~ PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
//~ PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
//~ PtrList<volScalarField> ghFluid(fluidRegions.size());
//~ PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
//~ PtrList<volScalarField> p_rghFluid(fluidRegions.size());
//~ PtrList<radiation::radiationModel> radiation(fluidRegions.size());

PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());

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
}
