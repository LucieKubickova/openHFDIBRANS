/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    scalarHFDIBTransportFoam

Group
    grpBasicSolvers

Description
    Passive scalar transport equation solver.

    \heading Solver details
    The equation is given by:

    \f[
        \ddt{T} + \div \left(\vec{U} T\right) - \div \left(D_T \grad T \right)
        = S_{T}
    \f]

    Where:
    \vartable
        T       | Passive scalar
        D_T     | Diffusion coefficient
        S_T     | Source
    \endvartable

    \heading Required fields
    \plaintable
        T       | Passive scalar
        U       | Velocity [m/s]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

#include "openHFDIBRANS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Passive scalar transport equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    // read simple dict
    dictionary HFDIBSIMPLEDict = simple.dict().subDict("HFDIB").subDict("T");
    word surfaceType;
    HFDIBSIMPLEDict.lookup("surfaceType") >> surfaceType;
    scalar boundaryVal = readScalar(HFDIBSIMPLEDict.lookup("boundaryValue"));
    scalar tolEqn = readScalar(HFDIBSIMPLEDict.lookup("tolEqn"));
    scalar maxEqnIters = readScalar(HFDIBSIMPLEDict.lookup("maxEqnIters"));
    scalar TIn = readScalar(HFDIBSIMPLEDict.lookup("valInside"));

    // prepare HFDIBRANS
    openHFDIBRANS HFDIBRANS(mesh, lambda);
    HFDIBRANS.createBaseSurface(surface, surfaceType, boundaryVal);

    // prepare T field
    Ti *= 0.0;
    Ti.correctBoundaryConditions();

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // HFDIBRANS -- NOTE: this before or inside non-orthogonal loop?
        HFDIBRANS.computeTi(T, Ti, surfaceType, TIn);
        HFDIBRANS.updateSurface(surface, surfaceType);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);

            for (label nCorr = 0; nCorr < maxEqnIters; nCorr++)
            {
                // HFDIBRANS: update source
                Tq = surface*(TEqn.A()*Ti - TEqn.H());
                solve(TEqn == Tq);

                Info << "HFDIB: Max error in T -> Ti is " << (max(mag(surface*(Ti - T))).value()) << endl;

                if (max(surface*(Ti - T)).value() < tolEqn)
                {
                    Info << "HFDIB: T converged to Ti within max tolerance " << tolEqn << endl;
                    break;
                }

                // apply correction
                T += 1.0*surface*(Ti - T);
            }

            fvOptions.correct(T);
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
