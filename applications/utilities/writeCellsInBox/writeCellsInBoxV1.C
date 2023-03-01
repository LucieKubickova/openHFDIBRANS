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
        "    xMin - minimal x coordinate\n"
        "    xMax - maximal x coordinate\n"
        "    yMin - minimal y coordinate\n"
        "    yMax - maximal y coordinate\n"
        "    zMin - minimal z coordinate\n"
        "    zMax - maximal z coordinate\n"
    );

    argList::noParallel();
    argList::validArgs.append("xMin");
    argList::validArgs.append("xMax");
    argList::validArgs.append("yMin");
    argList::validArgs.append("yMax");
    argList::validArgs.append("zMin");
    argList::validArgs.append("zMax");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // unpack arguments
    scalar xMin = -1*args.argRead<scalar>(1);
    scalar xMax = args.argRead<scalar>(2);
    scalar yMin = -1*args.argRead<scalar>(3);
    scalar yMax = args.argRead<scalar>(4);
    scalar zMin = -1*args.argRead<scalar>(5);
    scalar zMax = args.argRead<scalar>(6);

    word outDir("./constant/polyMesh/sets");
    if (!isDir(outDir))
    {
        mkDir(outDir);
    }

    OFstream outFile(outDir + "/cellsInBox");
    outFile << "cellI,x,y,z" << endl;
    
    // loop over all cells and write
    forAll(mesh.C(), cellI)
    {
        scalar x = mesh.C()[cellI].x();
        scalar y = mesh.C()[cellI].y();
        scalar z = mesh.C()[cellI].z();

        if (x >= xMin and x <= xMax and y >= yMin and y <= yMax and z >= zMin and z <= zMax)
        {
            outFile << cellI << "," << x << "," << y << "," << z << endl;
        }
    }

	Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
