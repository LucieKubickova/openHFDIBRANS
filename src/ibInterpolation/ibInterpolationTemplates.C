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

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::unifunctionalInterp
(
    ITstream& interpScheme,
    volTypeField& phi,
    volTypeField& phii,
    List<Type>& dirichletVals,
    List<scalar>& scales
)
{
    // create interpolator
    autoPtr<interpolation<Type>> interpPhi =
                   interpolation<Type>::New(HFDIBInterpDict_, phi);

    // read chosen interpolation function
    word interpType = interpScheme[1].wordToken();

    // prepare chosen interpolation function
    autoPtr<ibScheme> interpFunc = chosenInterpFunc(interpType);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // interpolate
        phii[cellI] = interpFunc->interpolate(phi, *interpPhi, body_, dirichletVals[bCell], scales[bCell], boundaryDists_[bCell].first(), intInfoList_[bCell], cellI);
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::switchedInterp
(
    ITstream& interpScheme,
    volTypeField& phi,
    volTypeField& phii,
    List<Type>& dirichletVals,
    List<scalar>& scales,
    volScalarField& yPlusi,
    scalar yPlusLam
)
{
    // create interpolator
    autoPtr<interpolation<Type>> interpPhi =
                   interpolation<Type>::New(HFDIBInterpDict_, phi);

    // read chosen interpolation functions
    word lowerType = interpScheme[1].wordToken();
    word higherType = interpScheme[2].wordToken();

    // prepare chosen interpolation functions
    autoPtr<ibScheme> lowerFunc = chosenInterpFunc(lowerType);
    autoPtr<ibScheme> higherFunc = chosenInterpFunc(higherType);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // switch based on yPlus value
        if (yPlusi[cellI] > yPlusLam)
        {
            phii[cellI] = higherFunc->interpolate(phi, *interpPhi, body_, dirichletVals[bCell], scales[bCell], boundaryDists_[bCell].first(), intInfoList_[bCell], cellI);
        }

        else
        {
            phii[cellI] = lowerFunc->interpolate(phi, *interpPhi, body_, dirichletVals[bCell], scales[bCell], boundaryDists_[bCell].first(), intInfoList_[bCell], cellI);
        }
    }
}

// ************************************************************************* //
