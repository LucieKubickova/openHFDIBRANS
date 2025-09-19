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
Type quadraticScheme::interpolateT
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
        return dirichletVal; // Note (LK): should call the constant interpolation, just dunno how to do it effectively
        //~ return constant<Type, volTypeField>(phi, interpPhi, dirichletVal, scale, bCell);
        //~ return linear<Type, volTypeField>(phi, interpPhi, dirichletVal, scale, bCell);
    }

    else if (order == 1)
    {
        // value in the interpolation point
        //~ Type phiP1 = interpPhi.interpolate(intInfo[1].iPoint_, intInfo[1].iCell_) - dirichletVal;
        Type phiP1 = phiPs[1] - dirichletVal;

        // distance between interpolation points
        scalar deltaR = mag(intInfo[1].iPoint_ - intInfo[0].iPoint_);

        // first polynomial coefficient
        Type linCoeff = phiP1/(deltaR+SMALL);

        // interpolated value
        return linCoeff*ds + dirichletVal;
    }

    // values in the interpolation points
    //~ Type phiP1 = interpPhi.interpolate(intInfo[1].iPoint_, intInfo[1].iCell_) - dirichletVal;
    //~ Type phiP2 = interpPhi.interpolate(intInfo[2].iPoint_, intInfo[2].iCell_) - dirichletVal;
    Type phiP1 = phiPs[1] - dirichletVal;
    Type phiP2 = phiPs[2] - dirichletVal;

    // distance between interpolation points
    scalar deltaR2 = mag(intInfo[2].iPoint_ - intInfo[1].iPoint_);
    scalar deltaR1 = mag(intInfo[1].iPoint_ - intInfo[0].iPoint_);

    // second polynomial coefficient
    Type quadCoeff = (phiP2 - phiP1)*deltaR1 - phiP1*deltaR2;
    quadCoeff      /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

    // first polynomial coefficient
    Type linCoeff  = (phiP1-phiP2)*Foam::pow(deltaR1,2.0);
    linCoeff  += 2.0*phiP1*deltaR1*deltaR2;
    linCoeff  += phiP1*Foam::pow(deltaR2,2.0);
    linCoeff  /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

    // interpolated value
    return quadCoeff*ds*ds + linCoeff*ds + dirichletVal;
}

// ************************************************************************* //
