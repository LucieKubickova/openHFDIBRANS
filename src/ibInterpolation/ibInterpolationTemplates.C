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

//~ #include "ibInterpolation.H"

#define ORDER 2

using namespace Foam;

//--------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
Type ibInterpolation::polynom
(
    volTypeField& phi,
    interpolation<Type>& interpPhi,
    Type& dirichletVal,
    label bCell
)
{
    // get boundary cell information
    label cellI = boundaryCells_[bCell].first();
    scalar ds = boundaryDists_[bCell].first();
    interpolationInfo intInfo = intInfoList_[bCell];

    // interpolated value
    Type toReturn;

    // interpolation
    switch(intInfo.order_)
    {
        case 0:
        {
            // interpolated value
            toReturn = body_[cellI]*dirichletVal + (1-body_[cellI])*phi[cellI];
            break;
        }

        case 1:
        {
            // value in the interpolation point
            Type phiP1 = interpPhi.interpolate(intInfo.intPoints_[1], intInfo.intCells_[0]) - dirichletVal;

            // distance between interpolation points
            scalar deltaR = mag(intInfo.intPoints_[1] - intInfo.intPoints_[0]);

            // first polynomial coefficient
            Type linCoeff = phiP1/(deltaR+SMALL);

            // interpolated value
            toReturn = linCoeff*ds + dirichletVal;
            break;
        }

        case 2:
        {
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
            toReturn = quadCoeff*ds*ds + linCoeff*ds + dirichletVal;
            break;
        }
    }

    // return
    return toReturn;
}

//--------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
Type ibInterpolation::logarithm
(
    volTypeField& phi,
    interpolation<Type>& interpPhi,
    scalar& dirichletVal,
    scalar& logScale,
    label bCell
)
{
    // get boundary cell information
    scalar ds = boundaryDists_[bCell].first();
    interpolationInfo intInfo = intInfoList_[bCell];

    // interpolated value
    Type toReturn;

    // interpolation
    switch(intInfo.order_)
    {
        case 0:
        {
            // interpolated value
            toReturn = polynom<Type, volTypeField>(phi, interpPhi, dirichletVal, bCell);
            break;
        }

        case 1: case 2: // what more can I do with 2 interpolation points?
        {
            // value in the first interpolation point
            Type phiP1 = interpPhi.interpolate(intInfo.intPoints_[1], intInfo.intCells_[0]) - dirichletVal;

            // distance between interpolation points
            scalar deltaR = mag(intInfo.intPoints_[1] - intInfo.intPoints_[0]);

            // compute A-log coefficient: 
            // y = A*ln(B*x + C) + D
            // y = A*ln(logScale*x + 1) + dirichletVal
            scalar B = logScale;
            scalar C = 1.0;
            Type A = phiP1/Foam::log(B*deltaR + C);

            // interpolated value
            toReturn = A*Foam::log(B*ds + C) + dirichletVal;
        }
    }

    // return
    return toReturn;
}

//--------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::polynomialInterp
(
    volTypeField& phi,
    volTypeField& phii,
    List<Type>& dirichletVals
)
{
    // create interpolator
    autoPtr<interpolation<Type>> interpPhi =
                   interpolation<Type>::New(HFDIBInterpDict_, phi);

    // interpolate and assign values to the imposed field based
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // interpolate and assign
        phii[cellI] = polynom<Type, volTypeField>(phi, *interpPhi, dirichletVals[bCell], bCell);
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::logarithmicInterp
(
    volTypeField& phi,
    volTypeField& phii,
    List<Type>& dirichletVals,
    List<scalar>& logScales
)
{
    // create interpolator
    autoPtr<interpolation<Type>> interpPhi =
                   interpolation<Type>::New(HFDIBInterpDict_, phi);

    // interpolate and assign values to the imposed field based
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // interpolate and assign
        phii[cellI] = logarithm<Type, volTypeField>(phi, *interpPhi, dirichletVals[bCell], logScales[bCell], bCell);
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::polySwitchLogInterp
(
    volTypeField& phi,
    volTypeField& phii,
    List<Type>& dirichletVals,
    List<scalar>& logScales,
    List<scalar>& yPlusi,
    scalar yPlusLam
)
{
    // create interpolator
    autoPtr<interpolation<Type>> interpPhi =
                   interpolation<Type>::New(HFDIBInterpDict_, phi);

    // interpolate and assign values to the imposed field based
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // switch base on yPlus value
        if (yPlusi[bCell] > yPlusLam)
        {
            // logarithmic interpolation
            phii[cellI] = logarithm<Type, volTypeField>(phi, *interpPhi, dirichletVals[bCell], logScales[bCell], bCell);
        }

        else
        {
            // polynomial interpolation
            phii[cellI] = polynom<Type, volTypeField>(phi, *interpPhi, dirichletVals[bCell], bCell);
        }
    }
}

// ************************************************************************* //
