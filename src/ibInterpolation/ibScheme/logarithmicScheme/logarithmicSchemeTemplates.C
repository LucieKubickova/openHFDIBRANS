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
        return body[cellI]*dirichletVal + (1-body[cellI])*phi[cellI]; // UGLYYYYYYYYYYYYYYYYYYY
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
