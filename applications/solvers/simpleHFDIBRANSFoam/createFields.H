Info << "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating field Ui\n" << endl;
volVectorField Ui
(
    IOobject
    (
        "Ui",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U
);

Info << "Reading field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating field surface\n" << endl;
volScalarField surface
(
    IOobject
    (
        "SIMPLE::surface",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,0)
);

Info << "Reading field f\n" << endl;
volVectorField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

pressureReference pressureReference(p, simple.dict());

mesh.schemes().setFluxRequired(p.name());

autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

autoPtr<incompressible::HFDIBmomentumTransportModel> turbulence
(
    incompressible::HFDIBmomentumTransportModel::New(U, phi, viscosity)
);

#include "readGravitationalAcceleration.H"
#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
