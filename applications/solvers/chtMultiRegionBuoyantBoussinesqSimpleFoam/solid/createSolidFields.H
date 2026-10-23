    // Initialise solid field pointer lists
    PtrList<volScalarField> TSolid(solidRegions.size());
    PtrList<volVectorField> USolid(solidRegions.size());
    PtrList<surfaceScalarField> phiSolid(solidRegions.size());
    PtrList<coordinateSystem> coordinates(solidRegions.size());
    PtrList<fv::options> solidHeatSources(solidRegions.size());
    PtrList<volScalarField> kappaSolid(solidRegions.size());
    PtrList<singlePhaseTransportModel> laminarTransportSolid(solidRegions.size());
    PtrList<uniformDimensionedScalarField> hRefSolid(solidRegions.size());
    PtrList<uniformDimensionedScalarField> TRefSolid(solidRegions.size());
    PtrList<uniformDimensionedScalarField> cpRefSolid(solidRegions.size());
    PtrList<uniformDimensionedScalarField> rhoRefSolid(solidRegions.size());
    PtrList<uniformDimensionedScalarField> betaSolid(solidRegions.size());
    PtrList<uniformDimensionedScalarField> alphaSolid(solidRegions.size());
    PtrList<radiation::radiationModel> radiations(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::options(solidRegions[i])
        );

        Info<< "    Adding to TSolid\n" << endl;
        TSolid.set  
        (           
            i,      
            new volScalarField
            (   
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );
    
        Info<< "    Adding to radiations\n" << endl;
        radiations.set(i, radiation::radiationModel::New(TSolid[i]));

        Info<< "    Adding to USolid\n" << endl;
        USolid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedVector("zero", dimLength/dimTime, vector::zero)
            )
        );      
                
        Info<< "    Adding to phiSolid\n" << endl;
        phiSolid.set
        (           
            i,      
            new surfaceScalarField
            (   
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(USolid[i])
                    & solidRegions[i].Sf()
            )
        );

        Info<< "    Adding to kappaSolid\n" << endl;
        kappaSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappa",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                //~ dimensionedScalar("kappa", dimLength*dimLength*dimLength*dimLength/dimTime/dimTime/dimTime/dimTemperature, 1.0)
                dimensionedScalar("kappa", dimMass*dimLength/dimTime/dimTime/dimTime/dimTemperature, 1.0)
            )
        );
    
        Info<< "    Adding to laminarTransportSolid\n" << endl;
        laminarTransportSolid.set
        (
            i,
            new singlePhaseTransportModel
            (
                USolid[i],
                phiSolid[i]
            )
        );
    
        Info<< "    Adding to TRefSolid\n" << endl;
        TRefSolid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "TRef",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar("TRef", dimTemperature, laminarTransportSolid[i])
            )
        );

        Info<< "    Adding to cpRefSolid\n" << endl;
        cpRefSolid.set
        (
            i, 
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "cpRef",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar("cpRef", dimLength*dimLength/dimTime/dimTime/dimTemperature, laminarTransportSolid[i])
            )
        );
    
        Info<< "    Adding to rhoRefSolid\n" << endl;
        rhoRefSolid.set
        (
            i, 
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "rhoRef",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar("rhoRef", dimMass/dimLength/dimLength/dimLength, laminarTransportSolid[i])
            )
        );
    
        Info<< "    Adding to betaSolid\n" << endl;
        betaSolid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "beta",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar("beta", dimless/dimTemperature, laminarTransportSolid[i])
            )
        );
    
        Info<< "    Adding to alphaSolid\n" << endl;
        alphaSolid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "alpha",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar("alpha", dimLength*dimLength/dimTime, laminarTransportSolid[i])
            )
        );
    }
