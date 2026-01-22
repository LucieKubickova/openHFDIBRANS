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
Type zeroGradientScheme::interpolateT
(
    volTypeField& phi,
    List<Type>& phiPs,
    Type& dirichletVal,
    scalar& scale,
    scalar& ds,
    List<intPoint>& intInfo,
    label& cellI
)
{
    // get interpolation order
    label order = getIntOrder(intInfo);

    // check whether there are enough interpolation points
    if (order == 0)
    {
        return 0.0*dirichletVal;
    }

    if (order == 1)
    {
        // value in the interpolation point
        //~ Type phiP1 = interpPhi.interpolate(intInfo[1].iPoint_, intInfo[1].iCell_);
        Type phiP1 = phiPs[1];

        // interpolated value
        return phiP1;
    }

    // values in the interpolation points
    //~ Type phiP1 = interpPhi.interpolate(intInfo[1].iPoint_, intInfo[1].iCell_);
    //~ Type phiP2 = interpPhi.interpolate(intInfo[2].iPoint_, intInfo[2].iCell_);
    Type phiP1 = phiPs[1];
    Type phiP2 = phiPs[2];

    // distance between interpolation points
    scalar deltaR2 = mag(intInfo[2].iPoint_ - intInfo[1].iPoint_);
    scalar deltaR1 = mag(intInfo[1].iPoint_ - intInfo[0].iPoint_);

    // second polynomial coefficient
    Type quadCoeff = phiP1 - phiP2;
    //~ quadCoeff      /= ((Foam::sqr(deltaR2) - Foam::sqr(deltaR1)) + SMALL);
    quadCoeff /= (-1*deltaR2*(deltaR2 + 2*deltaR1) + SMALL);

    // value at the delta = 0
    Type phiS = phiP1 - quadCoeff*Foam::sqr(deltaR1);

    // interpolated value
    return quadCoeff*Foam::sqr(ds) + phiS;
}

// ************************************************************************* //
