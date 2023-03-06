/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

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
        return body[cellI]*dirichletVal + (1-body[cellI])*phi[cellI]; // UGLYYYYYYYYYYYYYYYYYYYYY
        //~ return linear<Type, volTypeField>(phi, interpPhi, dirichletVal, scale, bCell);
    }

    else if (intInfo.order_ == 1)
    {
        // value in the interpolation point
        Type phiP1 = interpPhi.interpolate(intInfo.intPoints_[1], intInfo.intCells_[0]) - dirichletVal;

        // distance between interpolation points
        scalar deltaR = mag(intInfo.intPoints_[1] - intInfo.intPoints_[0]);

        // first polynomial coefficient
        Type linCoeff = phiP1/(deltaR+SMALL);

        // interpolated value
        return linCoeff*ds + dirichletVal;
        // UGLYYYYYYYYYYYYYYYYYYYY
    }

    // values in the interpolation points
    Type phiP1 = interpPhi.interpolate(intInfo.intPoints_[1], intInfo.intCells_[0]) - dirichletVal;
    Type phiP2 = interpPhi.interpolate(intInfo.intPoints_[2], intInfo.intCells_[1]) - dirichletVal;

    // distance between interpolation points
    scalar deltaR2 = mag(intInfo.intPoints_[2] - intInfo.intPoints_[1]);
    scalar deltaR1 = mag(intInfo.intPoints_[1] - intInfo.intPoints_[0]);

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
