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

    // read the k field
    Info << "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
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
    const tmp<volScalarField> tnut(turbulence->nut());
    const volScalarField& nut = tnut();

    // correct nut
    HFDIBRANS.updateUTau(k);
    HFDIBRANS.correctNut(const_cast<volScalarField&>(nut), k, const_cast<volScalarField&>(nu));

    // calculate wall shear stress
    HFDIBRANS.calculateWallShearStress(U, nu);

    // write fields
    volVectorField& tauwi = HFDIBRANS.getWallShearStress();
    tauwi.write();

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
