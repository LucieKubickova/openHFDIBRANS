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
#include "fvcSmooth.H"
#include "wallPolyPatch.H"
#include "faceSet.H"

int main(int argc, char *argv[])
{
    argList::addNote // help message
    (
        "Input arguments:\n"
        "----------------\n"
        "    surfVal - value at the surface\n"
    );

    argList::noParallel();
    argList::validArgs.append("surfVal");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

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

    // set parameters
    scalar thrSurf(readScalar(HFDIBDEMDict.lookup("surfaceThreshold")));

    // unpack arguments
    scalar surfVal = args.argRead<scalar>(1);

	// create the binary represented porous field and call it lambda
	volScalarField lambda
	(
		IOobject
		(
			"lambda",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("lambda",dimless,0.0)
	);

    // load lambda zone
    const word cellZoneName = "lambdaZone";
    const label zoneID = mesh.cellZones().findZoneID(cellZoneName);

    const labelList lambdaCells = mesh.cellZones()[zoneID];

    forAll(lambdaCells, lID)
    {
        lambda[lambdaCells[lID]] = 1.0;
    }

    forAll(lambdaCells, lID)
    {
        label cellID = lambdaCells[lID];

        forAll(mesh.cellCells()[cellID], nID)
        {
            if (lambda[mesh.cellCells()[cellID][nID]] == 0.0)
            {
                lambda[cellID] = surfVal;
            }
        }
    }

    // write the resulting fields
    lambda.write();

    // create a face set consisting boundary faces the adjacent cells of which have lambda > 0.5
    DynamicList<label> bFacesInsideLambda;

    forAll(mesh.boundaryMesh(), patchID)
    {
        if (isA<wallPolyPatch>(mesh.boundaryMesh()[patchID]))
        {
            forAll(mesh.boundary()[patchID].Cf(), faceID)
            {
                label owner = mesh.faceOwner()[mesh.boundary()[patchID].start() + faceID];

                if (lambda[owner] > thrSurf)
                {
                    bFacesInsideLambda.append(mesh.boundary()[patchID].start() + faceID);
                }
            }
        }
    }

    // write the face set
    faceSet inFacesSet
    (
        IOobject
        (
            "bFacesInsideLambda",
            polyMesh::meshSubDir/"sets",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    word outDir("./constant/polyMesh/sets");
    if (!isDir(outDir))
    {
        mkDir(outDir);
    }

    OFstream outFile(outDir + "/bFacesInsideLambda");
    inFacesSet.writeHeader(outFile);
    outFile << bFacesInsideLambda << endl;
    outFile << "// ************************************************************************* //";

	Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
