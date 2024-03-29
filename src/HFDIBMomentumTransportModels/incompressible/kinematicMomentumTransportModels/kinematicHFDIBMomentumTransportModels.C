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

#include "kinematicHFDIBMomentumTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeBaseHFDIBMomentumTransportModel
(
    geometricOneField,
    geometricOneField,
    incompressibleHFDIBMomentumTransportModel,
    IncompressibleHFDIBMomentumTransportModel,
    transportModel
);


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "../../momentumTransportModels/lnInclude/HFDIBStokes.H"
makeHFDIBLaminarModel(HFDIBStokes);

#include "../../momentumTransportModels/lnInclude/HFDIBGeneralizedNewtonian.H"
makeHFDIBLaminarModel(HFDIBGeneralizedNewtonian);

#include "../../momentumTransportModels/lnInclude/HFDIBMaxwell.H"
makeHFDIBLaminarModel(HFDIBMaxwell);

#include "../../momentumTransportModels/lnInclude/HFDIBGiesekus.H"
makeHFDIBLaminarModel(HFDIBGiesekus);

#include "../../momentumTransportModels/lnInclude/HFDIBPTT.H"
makeHFDIBLaminarModel(HFDIBPTT);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "../../momentumTransportModels/lnInclude/HFDIBKOmega.H"
makeHFDIBRASModel(HFDIBKOmega);

#include "../../momentumTransportModels/lnInclude/HFDIBKOmegaSST.H"
makeHFDIBRASModel(HFDIBKOmegaSST);

#include "../../momentumTransportModels/lnInclude/HFDIBKEpsilon.H"
makeHFDIBRASModel(HFDIBKEpsilon);

#include "../../momentumTransportModels/lnInclude/HFDIBRealizableKE.H"
makeHFDIBRASModel(HFDIBRealizableKE);

// ************************************************************************* //
