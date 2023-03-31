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
openHFDIBRANS is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIBRANS. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with RAS turbulence modeling approach and
    wall functions (Kubickova and Isoz, 2023)

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicHFDIBMomentumTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "addToRunTimeSelectionTable.H"
#include "HFDIBMomentumTransportModel.H"
#include "IncompressibleHFDIBMomentumTransportModel.H"
#include "transportModel.H"
#include "HFDIBRASModel.H"

#include "triSurfaceMesh.H"

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

    Info << "\nStarting time loop\n" << endl;

    // read simple dict
    dictionary HFDIBSIMPLEDict = simple.dict().subDict("HFDIB");
    word surfaceType;
    HFDIBSIMPLEDict.lookup("surfaceType") >> surfaceType;
    scalar boundaryVal = readScalar(HFDIBSIMPLEDict.lookup("boundaryValue"));
    scalar tolUEqn = readScalar(HFDIBSIMPLEDict.lookup("tolUEqn"));
    scalar maxUEqnIters = readScalar(HFDIBSIMPLEDict.lookup("maxUEqnIters"));

    // prepare HFDIBRANS
    openHFDIBRANS HFDIBRANS(mesh, lambda, simulationType);
    HFDIBRANS.createBaseSurface(surface, surfaceType, boundaryVal);

    Ui *= 0.0;
    Ui.correctBoundaryConditions();

    while (simple.loop(runTime))
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "pEqn.H"

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
