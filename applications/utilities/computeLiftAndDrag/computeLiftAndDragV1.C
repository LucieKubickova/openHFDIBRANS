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
#include "cellSet.H"
#include "unitConversion.H"

int main(int argc, char *argv[])
{
    argList::addNote // help message
    (
        "Input arguments:\n"
        "----------------\n"
        "    alpha - angle of attack in degrees\n"
        "    cellSetName - name of the cell set to work on\n"
    );

    argList::noParallel();
    argList::validArgs.append("alpha");
    argList::validArgs.append("cellSetName");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // unpack arguments
    scalar alpha = degToRad(args.argRead<scalar>(1));
    word cellSetName = args.argRead<word>(2);

    // read the force field field
    Info << "Reading field f\n" << endl;
    volVectorField f
    (
        IOobject
        (
            "f",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // load the cell set
    cellSet boundaryCells
    (
        IOobject
        (
            cellSetName,
            "constant/polyMesh/sets",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // compute force
    vector F;
    scalar V(0.0);

    forAll(mesh.C(), cellI)
    {
        if (boundaryCells[cellI])
        {
            F += f[cellI]*mesh.V()[cellI];
            V += mesh.V()[cellI];
        }
    }

    // scale the force
    F /= V;

    // compute lift 2D
    scalar L = -1*F.x()*Foam::sin(alpha) + F.y()*Foam::cos(alpha);

    // compute drag 2D
    scalar D = F.x()*Foam::cos(alpha) + F.y()*Foam::sin(alpha);

    // print
    Info << "The resulting forces are" << endl;
    Info << "\tTotal: F = " << F << endl;
    Info << "\t Lift: L = " << L << endl;
    Info << "\t Drag: D = " << D << endl;
	Info << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
