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
                   interpolation<Type>::New(HFDIBInnerSchemes_, phi);

    // read chosen interpolation function
    word ibSchemeType = ibSchemeName[1].wordToken();

    // prepare chosen interpolation function
    autoPtr<ibScheme> interpFunc = chosenInterpFunc(ibSchemeType);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // interpolate
        phii[cellI] = interpFunc->interpolate
        (
            phi,
            *interpPhi,
            body_,
            dirichletVals[bCell],
            scales[bCell],
            boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_,
            lineIntInfoBoundary_->getIntPoints()[bCell],
            cellI
        );
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::lambdaBasedInterp
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
                   interpolation<Type>::New(HFDIBInnerSchemes_, phi);

    // read chosen interpolation function
    word ibSchemeType = ibSchemeName[1].wordToken();

    // prepare chosen interpolation function
    autoPtr<ibScheme> interpFunc = chosenInterpFunc(ibSchemeType);

    // interpolate and assign values to the imposed field
    forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    {
        // get cell label
        label cellI = surfaceCells_[Pstream::myProcNo()][sCell];

        // interpolate
        phii[cellI] = interpFunc->interpolate
        (
            phi,
            *interpPhi,
            body_,
            dirichletVals[sCell],
            scales[sCell],
            surfaceDists_[Pstream::myProcNo()][sCell],
            lineIntInfoSurface_->getIntPoints()[sCell],
            cellI
        );
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
                   interpolation<Type>::New(HFDIBInnerSchemes_, phi);

    // read chosen interpolation functions
    word lowerType = ibSchemeName[1].wordToken();
    word higherType = ibSchemeName[2].wordToken();

    // prepare chosen interpolation functions
    autoPtr<ibScheme> lowerFunc = chosenInterpFunc(lowerType);
    autoPtr<ibScheme> higherFunc = chosenInterpFunc(higherType);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // switch based on yPlus value
        if (yPlusi[cellI] > yPlusLam)
        {
            phii[cellI] = higherFunc->interpolate
            (
                phi,
                *interpPhi,
                body_,
                dirichletVals[bCell],
                scales[bCell],
                boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_,
                lineIntInfoBoundary_->getIntPoints()[bCell],
                cellI
            );
        }

        else
        {
            phii[cellI] = lowerFunc->interpolate
            (
                phi,
                *interpPhi,
                body_,
                dirichletVals[bCell],
                scales[bCell],
                boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_,
                lineIntInfoBoundary_->getIntPoints()[bCell],
                cellI
            );
        }
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::outerInnerInterp
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
                   interpolation<Type>::New(HFDIBInnerSchemes_, phi);

    // read chosen interpolation functions
    word outerType = ibSchemeName[1].wordToken();
    word innerType = ibSchemeName[2].wordToken();

    // prepare chosen interpolation functions
    autoPtr<ibScheme> outerFunc = chosenInterpFunc(outerType);
    autoPtr<ibScheme> innerFunc = chosenInterpFunc(innerType);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get outer cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // switch based on yPlus value
        if (yPlusi[outCellI] > yPlusLam)
        {
            // get the inner cell label
            label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

            // create new interpolation info to interpolate inside
            //~ interpolationInfo intInfoToSend(intInfoListBoundary_[Pstream::myProcNo()][bCell]);
            //~ intInfoToSend.intPoints_[1] = mesh_.C()[outCellI];
            //~ intInfoToSend.intPoints_[2] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intPoints_[1];
            //~ intInfoToSend.intCells_[0] = outCellI;
            //~ intInfoToSend.intCells_[1] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intCells_[0];

            // NOTE: logarithm of negative number?
            phii[inCellI] = innerFunc->interpolate
            (
                phi,
                *interpPhi,
                body_,
                dirichletVals[bCell],
                scales[bCell],
                boundaryCells_[Pstream::myProcNo()][bCell].sigma_,
                lineIntInfoBoundary_->getIntPoints()[bCell], // NOTE: wrong interpolation info passed
                outCellI
            );
            // NOTE: constant won't work
        }

        else
        {
            phii[outCellI] = outerFunc->interpolate
            (
                phi,
                *interpPhi,
                body_,
                dirichletVals[bCell],
                scales[bCell],
                boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_,
                lineIntInfoBoundary_->getIntPoints()[bCell],
                outCellI
            );
        }
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::innerInterp
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
                   interpolation<Type>::New(HFDIBInnerSchemes_, phi);

    // read chosen interpolation function
    word ibSchemeType = ibSchemeName[1].wordToken();

    // prepare chosen interpolation function
    autoPtr<ibScheme> interpFunc = chosenInterpFunc(ibSchemeType);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the inner cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

        // create new interpolation info to interpolate inside
        //~ interpolationInfo intInfoToSend(intInfoListBoundary_[Pstream::myProcNo()][bCell]);
        //~ intInfoToSend.intPoints_[1] = mesh_.C()[outCellI];
        //~ intInfoToSend.intPoints_[2] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intPoints_[1];
        //~ intInfoToSend.intCells_[0] = outCellI;
        //~ intInfoToSend.intCells_[1] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intCells_[0];
        
        // NOTE: logarithm of negative number?
        phii[inCellI] = interpFunc->interpolate
        (
            phi,
            *interpPhi,
            body_,
            dirichletVals[bCell],
            scales[bCell],
            boundaryCells_[Pstream::myProcNo()][bCell].sigma_,
            lineIntInfoBoundary_->getIntPoints()[bCell], // NOTE: wrong interpolation info passed
            outCellI
        );
    }
}

// ************************************************************************* //
