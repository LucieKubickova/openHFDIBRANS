/*---------------------------------------------------------------------------*\
                        _   _ ____ ____ _____ _____ _____ _____ _    _  _____
                       | | | |  __|  _ \_   _|  __ \  __ \  _  \ \  | |/  _  \
  ___  _ __   ___ _ __ | |_| | |_ | | | || | | |_/ / |_/ / |_| |  \ | |  |_|_/
 / _ \| '_ \ / _ \ '_ \|  _  |  _|| | | || | |  __ \  _ ||  _  | \ \| |\___  \
| (_) | |_) |  __/ | | | | | | |  | |/ / | |_| |_/ / | \ \ | | | |\ \ |/ |_|  |
 \___/| .__/ \___|_| |_\_| |_\_|  |___/ \___/\____/|_/  \_|| |_|_| \__|\_____/
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                    with R eynolds A veraged N avier S tokes equations          
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
    buoyantBoussinesqSimpleFoam

Group
    grpHeatTransferSolvers

Description
    Steady-state solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "HFDIBTurbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "HFDIBTurbulenceModel.H"
#include "IncompressibleHFDIBTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "HFDIBRASModel.H"
#include "HFDIBLaminarModel.H"

#include "triSurfaceMesh.H"
#include "openHFDIBRANS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for buoyant, turbulent flow"
        " of incompressible fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    // read simple dict for U
    dictionary HFDIBSIMPLEDictU = simple.dict().subDict("HFDIB").subDict("U");
    word surfaceTypeU;
    HFDIBSIMPLEDictU.lookup("surfaceType") >> surfaceTypeU;
    scalar boundaryValU = readScalar(HFDIBSIMPLEDictU.lookup("boundaryValue"));
    scalar tolEqnU = readScalar(HFDIBSIMPLEDictU.lookup("tolEqn"));
    scalar maxEqnItersU = readScalar(HFDIBSIMPLEDictU.lookup("maxEqnIters"));

    // read dict for T
    dictionary HFDIBSIMPLEDictT = simple.dict().subDict("HFDIB").subDict("T");
    word surfaceTypeT;
    HFDIBSIMPLEDictT.lookup("surfaceType") >> surfaceTypeT;
    scalar boundaryValT = readScalar(HFDIBSIMPLEDictT.lookup("boundaryValue"));
    scalar tolEqnT = readScalar(HFDIBSIMPLEDictT.lookup("tolEqn"));
    scalar maxEqnItersT = readScalar(HFDIBSIMPLEDictT.lookup("maxEqnIters"));
    scalar TIn = readScalar(HFDIBSIMPLEDictT.lookup("valInside"));

    // prepare HFDIBRANS
    openHFDIBRANS HFDIBRANS(mesh, lambda);
    HFDIBRANS.createBaseSurface(surfaceU, surfaceTypeU, boundaryValU);
    HFDIBRANS.createBaseSurface(surfaceT, surfaceTypeT, boundaryValT);

    Ui *= 0.0;
    Ui.correctBoundaryConditions();
    Ti *= 0.0;
    Ti.correctBoundaryConditions();

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "TEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct(HFDIBRANS);

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
