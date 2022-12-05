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
#include "OFstream.H"
#include "RASModel.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // read the omega field
    Info << "Reading field omega\n" << endl;
    volScalarField omega
    (
        IOobject
        (
            "omega",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // NOT USED NOW
    // read the epsilon field
    //~ Info << "Reading field epsilon\n" << endl;
    //~ volScalarField epsilon
    //~ (
        //~ IOobject
        //~ (
            //~ "epsilon",
            //~ runTime.timeName(),
            //~ mesh,
            //~ IOobject::READ_IF_PRESENT,
            //~ IOobject::NO_WRITE
        //~ ),
        //~ mesh
    //~ );

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

    // read the nut field
    Info << "Reading field nut\n" << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // read the nut field
    Info << "Reading field G\n" << endl;
    volScalarField::Internal G
    (
        IOobject
        (
            "G",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // prepare space to save
    word saveDir = mesh.time().rootPath() + "/" + mesh.time().globalCaseName() + "/postProcessing/readAndWriteTurbulenceFields";
    if(!isDir(saveDir))
        mkDir(saveDir);

    // prepare save file
    word fileName = "turbulenceFields.dat";

    autoPtr<OFstream> outFilePtr;
    outFilePtr.reset(new OFstream(saveDir/fileName));
    outFilePtr() << "x,y,z,omega,G,k,nut" << endl;

    // save data
    forAll(mesh.C(), cellI)
    {
        scalar x = mesh.C()[cellI].x();
        scalar y = mesh.C()[cellI].y();
        scalar z = mesh.C()[cellI].z();
        scalar omegai = omega[cellI];
        scalar Gi = G[cellI];
        scalar ki = k[cellI];
        scalar nuti = nut[cellI];

        outFilePtr() << x << "," << y << "," << z << "," << omegai << "," << Gi << "," << ki << "," << nuti << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
