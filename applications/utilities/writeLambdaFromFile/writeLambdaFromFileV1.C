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
        "    fileName - file with cell IDs and distance\n"
    );

    argList::noParallel();
    argList::validArgs.append("fileName");

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
    scalar intSpan(readScalar(HFDIBDEMDict.lookup("interfaceSpan")));

    // unpack arguments
    word fileName = args.argRead<word>(1);

    // compute average cell volume
    scalar VAve = 0.0;
    forAll(mesh.V(), i)
    {
        VAve += mesh.V()[i];
    }
    VAve /= mesh.V().size();

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

    // read data
    IOdictionary distanceDict
    (
        IOobject
        (
            fileName,
            runTime.constant() + "/polyMesh/sets",
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    List<label> cellList;
    List<scalar> distList;

    distanceDict.lookup("nearestCells") >> cellList;
    distanceDict.lookup("distances") >> distList;

    // set lambda
    forAll(cellList, i)
    {
        // get the cell label
        label cellI = cellList[i];

        // get the distance
        scalar ds = distList[i];

        // compute lambda
        lambda[cellI] = (1 - Foam::tanh(ds*intSpan/Foam::pow(VAve,0.333)))*0.5;

        if (lambda[cellI] < thrSurf)
        {
            lambda[cellI] = 0.0;
        }

        else if (lambda[cellI] > (1-thrSurf))
        {
            lambda[cellI] = 1.0;
        }
    }

    // write the resulting fields
    lambda.write();

	Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
