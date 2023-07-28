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

//~ #include "ibInterpolation.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::unifunctionalInterp
(
    ITstream& ibSchemeName,
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
    word ibSchemeType = ibSchemeName[1].wordToken();

    // prepare chosen interpolation function
    autoPtr<ibScheme> interpFunc = chosenInterpFunc(ibSchemeType);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // interpolate
        if (boundaryDists_[bCell].first() > 0.0)
        {
            phii[cellI] = interpFunc->interpolate(phi, *interpPhi, body_, dirichletVals[bCell], scales[bCell], boundaryDists_[bCell].first(), intInfoList_[bCell], cellI);
        }

        else
        {
            phii[cellI] = dirichletVals[bCell];
        }
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::switchedInterp
(
    ITstream& ibSchemeName,
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
    word lowerType = ibSchemeName[1].wordToken();
    word higherType = ibSchemeName[2].wordToken();

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
            if (boundaryDists_[bCell].first() > 0.0)
            {
                phii[cellI] = higherFunc->interpolate(phi, *interpPhi, body_, dirichletVals[bCell], scales[bCell], boundaryDists_[bCell].first(), intInfoList_[bCell], cellI);
            }

            else
            {
                phii[cellI] = dirichletVals[bCell];
            }
        }

        else
        {
            if (boundaryDists_[bCell].first() > 0.0)
            {
                phii[cellI] = lowerFunc->interpolate(phi, *interpPhi, body_, dirichletVals[bCell], scales[bCell], boundaryDists_[bCell].first(), intInfoList_[bCell], cellI);
            }

            else
            {
                phii[cellI] = dirichletVals[bCell];
            }
        }
    }
}
// ************************************************************************* //
