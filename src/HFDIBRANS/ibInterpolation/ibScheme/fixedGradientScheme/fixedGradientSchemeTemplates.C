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
Type fixedGradientScheme::interpolateT
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
    //~ if (intInfo.order_ == 0)
    //~ {
        //~ return dirichletVal; // UGLYYYYYYYYYYYYYYYYYYYYY
        //~ return linear<Type, volTypeField>(phi, interpPhi, dirichletVal, scale, bCell);
    //~ }

    if (intInfo.order_ == 1)
    {
        // value in the interpolation point
        Type phiP1 = interpPhi.interpolate(intInfo.intPoints_[1], intInfo.intCells_[0]);

        // distance between interpolation points
        scalar deltaR = mag(intInfo.intPoints_[1] - intInfo.intPoints_[0]);

        // value at the delta = 0
        Type phiS = phiP1 - dirichletVal*deltaR;

        // interpolated value
        return dirichletVal*deltaR + phiS;
    }

    // values in the interpolation points
    Type phiP1 = interpPhi.interpolate(intInfo.intPoints_[1], intInfo.intCells_[0]);
    Type phiP2 = interpPhi.interpolate(intInfo.intPoints_[2], intInfo.intCells_[1]);

    // distance between interpolation points
    scalar deltaR2 = mag(intInfo.intPoints_[2] - intInfo.intPoints_[1]);
    scalar deltaR1 = mag(intInfo.intPoints_[1] - intInfo.intPoints_[0]);

    // second polynomial coefficient
    Type quadCoeff = phiP1 - phiP2 + dirichletVal*(deltaR2);
    //~ quadCoeff      /= ((Foam::sqr(deltaR2) - Foam::sqr(deltaR1)) + SMALL);
    quadCoeff /= (-1*deltaR2*(deltaR2 + 2*deltaR1) + SMALL);

    // value at the delta = 0
    Type phiS = phiP1 - quadCoeff*Foam::sqr(deltaR1) - dirichletVal*deltaR1;

    // interpolated value
    return quadCoeff*Foam::sqr(ds) + dirichletVal*ds + phiS;
}

// ************************************************************************* //
