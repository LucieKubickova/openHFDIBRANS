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

using namespace Foam;

//--------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
Type logarithmicScheme::interpolateT
(
    volTypeField& phi,
    interpolation<Type>& interpPhi,
    const volScalarField& body,
    Type& dirichletVal,
    scalar& scale,
    scalar& ds,
    interpolationInfo& intInfo,
    label& cellI
)
{
    // check whether there are enough interpolation points
    if (intInfo.order_ == 0)
    {
        return dirichletVal; // UGLYYYYYYYYYYYYYYYYYYY
        //~ return constant<Type, volTypeField>(phi, interpPhi, dirichletVal, scale, bCell);
    }

    // value in the first interpolation point
    Type phiP1 = interpPhi.interpolate(intInfo.intPoints_[1], intInfo.intCells_[0]) - dirichletVal;

    // distance between interpolation points
    scalar deltaR = mag(intInfo.intPoints_[1] - intInfo.intPoints_[0]);

    // compute A-log coefficient: 
    // y = A*ln(B*x + C) + D
    // y = A*ln(logScale*x + 1) + dirichletVal
    scalar B = scale;
    scalar C = 1.0;
    Type A = phiP1/Foam::log(B*deltaR + C);

    // interpolated value
    return A*Foam::log(B*ds + C) + dirichletVal;
}

// ************************************************************************* //
