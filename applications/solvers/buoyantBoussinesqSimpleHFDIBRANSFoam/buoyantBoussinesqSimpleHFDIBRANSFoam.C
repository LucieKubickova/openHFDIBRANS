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
    openHFDIBRANS is licensed under the GNU LESSER GENERAL PUBLIC LICENSE
    (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the
    terms and conditions of version 3 of the GNU General Public License,
    supplemented by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIBRANS. If not, see
    <http://www.gnu.org/licenses/lgpl.html>.

Application
    buoyantBoussinesqSimpleFoam

Group
    grpHeatTransferSolvers

Description
    Steady-state solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho \f$ = the effective (driving) density
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
#include "kinematicHFDIBMomentumTransportModel.H"
#include "radiationModel.H"
#include "noRadiation.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "HFDIBMomentumTransportModel.H"
#include "IncompressibleHFDIBMomentumTransportModel.H"
#include "transportModel.H"
#include "HFDIBRASModel.H"
#include "HFDIBLaminarModel.H"

#include "triSurfaceMesh.H"
#include "openHFDIBRANS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // read simple dict for U
    dictionary HFDIBSIMPLEDictU = simple.dict().subDict("HFDIB").subDict("U");
    word surfaceTypeU;
    HFDIBSIMPLEDictU.lookup("surfaceType") >> surfaceTypeU;
    scalar boundaryValU = readScalar(HFDIBSIMPLEDictU.lookup("boundaryValue"));
    vector UIn = HFDIBSIMPLEDictU.lookupOrDefault<vector>("valInside", vector::zero);

    // read dict for T
    dictionary HFDIBSIMPLEDictT = simple.dict().subDict("HFDIB").subDict("T");
    word surfaceTypeT;
    HFDIBSIMPLEDictT.lookup("surfaceType") >> surfaceTypeT;
    scalar boundaryValT = readScalar(HFDIBSIMPLEDictT.lookup("boundaryValue"));
    scalar TIn = readScalar(HFDIBSIMPLEDictT.lookup("valInside"));

    // prepare HFDIBRANS
    openHFDIBRANS HFDIBRANS(mesh, lambda);
    HFDIBRANS.createBaseSurface(surfaceU, surfaceTypeU, boundaryValU);
    surfaceU.correctBoundaryConditions();
    HFDIBRANS.createBaseSurface(surfaceT, surfaceTypeT, boundaryValT);
    surfaceT.correctBoundaryConditions();

    Ui *= 0.0;
    Ui.correctBoundaryConditions();
    Ti *= 0.0;
    Ti.correctBoundaryConditions();

    Info << "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
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

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
