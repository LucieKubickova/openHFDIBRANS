/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicHFDIBMomentumTransportModel.H"
#include "fvOptions.H"

#include "addToRunTimeSelectionTable.H"
#include "HFDIBMomentumTransportModel.H"
#include "IncompressibleHFDIBMomentumTransportModel.H"
#include "transportModel.H"
#include "HFDIBRASModel.H"
#include "HFDIBLaminarModel.H"

#include "triSurfaceMesh.H"
#include "openHFDIBRANS.H"

#include "OFstream.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // read the U field
    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // read the p field
    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // read the phi field
    Info<< "Reading field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // read the lambda field
    Info << "Reading field lambda\n" << endl;
    volScalarField lambda
    (
        IOobject
        (
            "lambda",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // prepare fN field
    volVectorField fN
    (
        IOobject
        (
            "fN",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    );

    // prepare tauw field
    volVectorField tauw
    (
        IOobject
        (
            "tauw",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    );

    // prepare fT field
    volVectorField fT
    (
        IOobject
        (
            "fT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    );

    // load HFDIBDEM dictionary
    IOdictionary HFDIBDEMDict
    (
        IOobject
        (
            "HFDIBDEMDict",
            "constant",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // prepare turbulence model
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::HFDIBmomentumTransportModel> turbulence
    (
        incompressible::HFDIBmomentumTransportModel::New(U, phi, laminarTransport)
    );

    // prepare HFDIBRANS
    openHFDIBRANS HFDIBRANS(mesh, lambda);

    // read viscosity fields
    const tmp<volScalarField> tnu(turbulence->nu());
    const volScalarField& nu = tnu();

    // prepare dummy fields
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength/dimTime, 0.0)
    );

    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength/dimTime/dimTime, 0.0)
    );

    // correct nut
    HFDIBRANS.correctNut(nut, k, const_cast<volScalarField&>(nu));

    // calculate wall shear stress
    HFDIBRANS.calculateWallShearStress(tauw, U, nu);

    // calculate surface area
    HFDIBRANS.calculateSurfaceArea();

    // prepare dictionary
    dictionary forceDict = HFDIBDEMDict.subDict("forceCoeffs");

    // calculate forces
    HFDIBRANS.calculateForces(fN, fT, tauw, p, forceDict);

    // calculate lift and drag
    scalar Cl(0.0);
    scalar Cd(0.0);
    HFDIBRANS.calculateForceCoeffs(Cl, Cd, fN, fT, forceDict);

    // write
    Info << "Force coeffs are:" << endl;
    Info << "    Cl = " << Cl << endl;
    Info << "    Cd = " << Cd << endl;

    Info << endl << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
