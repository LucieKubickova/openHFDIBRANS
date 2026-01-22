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

    // read values in interpolation points 
    List<List<intPoint>>& intInfos = lineIntInfoBoundary_->getIntPoints();

    // get values in interpolation points
    List<List<Type>> phiInIntPoints(intInfos.size());
    interpolateToIntPoints<Type, volTypeField>(phi, *interpPhi, intInfos, phiInIntPoints);

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // interpolate
        phii[cellI] = interpFunc->interpolate
        (
            phi,
            phiInIntPoints[bCell],
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

    // read values in interpolation points 
    // Note (LK): surface int info not checked in parallel
    List<List<intPoint>>& intInfos = lineIntInfoSurface_->getIntPoints();

    // get values in interpolation points
    List<List<Type>> phiInIntPoints(intInfos.size());
    interpolateToIntPoints<Type, volTypeField>(phi, *interpPhi, intInfos, phiInIntPoints);

    // interpolate and assign values to the imposed field
    forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    {
        // get cell label
        label cellI = surfaceCells_[Pstream::myProcNo()][sCell].sCell_;

        // interpolate
        phii[cellI] = interpFunc->interpolate
        (
            phi,
            phiInIntPoints[sCell],
            dirichletVals[sCell],
            scales[sCell],
            surfaceCells_[Pstream::myProcNo()][sCell].sigma_,
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

    // read values in interpolation points 
    List<List<intPoint>>& intInfos = lineIntInfoBoundary_->getIntPoints();

    // get values in interpolation points
    List<List<Type>> phiInIntPoints(intInfos.size());
    interpolateToIntPoints<Type, volTypeField>(phi, *interpPhi, intInfos, phiInIntPoints);

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
                phiInIntPoints[bCell],
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
                phiInIntPoints[bCell],
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

    // prepare sync
    List<DynamicList<label>> labelsToSync(Pstream::nProcs());
    List<DynamicList<Type>> phisToSync(Pstream::nProcs());

    // read values in interpolation points 
    List<List<intPoint>>& intInfos = lineIntInfoBoundary_->getIntPoints();

    // get values in interpolation points
    List<List<Type>> phiInIntPoints(intInfos.size());
    interpolateToIntPoints<Type, volTypeField>(phi, *interpPhi, intInfos, phiInIntPoints);

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
            label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;

            // create new interpolation info to interpolate inside
            //~ interpolationInfo intInfoToSend(intInfoListBoundary_[Pstream::myProcNo()][bCell]);
            //~ intInfoToSend.intPoints_[1] = mesh_.C()[outCellI];
            //~ intInfoToSend.intPoints_[2] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intPoints_[1];
            //~ intInfoToSend.intCells_[0] = outCellI;
            //~ intInfoToSend.intCells_[1] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intCells_[0];

            // Note (LK): logarithm of negative number?
            // Note (LK): constant won't work
            Type phiS = innerFunc->interpolate
            (
                phi,
                phiInIntPoints[bCell],
                dirichletVals[bCell],
                scales[bCell],
                boundaryCells_[Pstream::myProcNo()][bCell].sigma_,
                lineIntInfoBoundary_->getIntPoints()[bCell], // Note (LK): wrong interpolation info passed, but not used now
                outCellI
            );

            if (Pstream::myProcNo() == iProc)
            {
                phii[inCellI] = phiS;
            }
            else
            {
                // save to sync
                labelsToSync[iProc].append(inCellI);
                phisToSync[iProc].append(phiS);
            }
        }

        else
        {
            phii[outCellI] = outerFunc->interpolate
            (
                phi,
                phiInIntPoints[bCell],
                dirichletVals[bCell],
                scales[bCell],
                boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_,
                lineIntInfoBoundary_->getIntPoints()[bCell],
                outCellI
            );
        }
    }

    // sync with other processors
    PstreamBuffers pBufsLabels(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsPhis(Pstream::commsTypes::nonBlocking);

    // send
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendLabels(proci, pBufsLabels);
            UOPstream sendPhis(proci, pBufsPhis);
            sendLabels << labelsToSync[proci];
            sendPhis << phisToSync[proci];
        }
    }
    
    pBufsLabels.finishedSends();
    pBufsPhis.finishedSends();
    
    // recieve
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvLabels(proci, pBufsLabels);
            UIPstream recvPhis(proci, pBufsPhis);
            DynamicList<label> recLabels (recvLabels);
            DynamicList<Type> recPhis (recvPhis);

            forAll(recLabels, rCell)
            {
                // get the cell label
                label cellI = recLabels[rCell];

                // assign phi
                phii[cellI] = recPhis[rCell];
            }
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

    // read values in interpolation points 
    List<List<intPoint>>& intInfos = lineIntInfoBoundary_->getIntPoints();

    // get values in interpolation points
    List<List<Type>> phiInIntPoints(intInfos.size());
    interpolateToIntPoints<Type, volTypeField>(phi, *interpPhi, intInfos, phiInIntPoints);

    // prepare sync
    List<DynamicList<label>> labelsToSync(Pstream::nProcs());
    List<DynamicList<Type>> phisToSync(Pstream::nProcs());

    // interpolate and assign values to the imposed field
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the inner cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;
        label iProc = boundaryCells_[Pstream::myProcNo()][bCell].iProc_;

        // create new interpolation info to interpolate inside
        //~ interpolationInfo intInfoToSend(intInfoListBoundary_[Pstream::myProcNo()][bCell]);
        //~ intInfoToSend.intPoints_[1] = mesh_.C()[outCellI];
        //~ intInfoToSend.intPoints_[2] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intPoints_[1];
        //~ intInfoToSend.intCells_[0] = outCellI;
        //~ intInfoToSend.intCells_[1] = intInfoListBoundary_[Pstream::myProcNo()][bCell].intCells_[0];
        
        // Note (LK): logarithm of negative number?
        Type phiS = interpFunc->interpolate
        (
            phi,
            phiInIntPoints[bCell],
            dirichletVals[bCell],
            scales[bCell],
            boundaryCells_[Pstream::myProcNo()][bCell].sigma_,
            lineIntInfoBoundary_->getIntPoints()[bCell], // Note (LK): wrong interpolation info passed, but not used now
            outCellI
        );

        if (Pstream::myProcNo() == iProc)
        {
            phii[inCellI] = phiS;
        }
        else
        {
            // save to sync
            labelsToSync[iProc].append(inCellI);
            phisToSync[iProc].append(phiS);
        }
    }

    // sync with other processors
    PstreamBuffers pBufsLabels(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsPhis(Pstream::commsTypes::nonBlocking);

    // send
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendLabels(proci, pBufsLabels);
            UOPstream sendPhis(proci, pBufsPhis);
            sendLabels << labelsToSync[proci];
            sendPhis << phisToSync[proci];
        }
    }
    
    pBufsLabels.finishedSends();
    pBufsPhis.finishedSends();
    
    // recieve
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvLabels(proci, pBufsLabels);
            UIPstream recvPhis(proci, pBufsPhis);
            DynamicList<label> recLabels (recvLabels);
            DynamicList<Type> recPhis (recvPhis);

            forAll(recLabels, rCell)
            {
                // get the cell label
                label cellI = recLabels[rCell];

                // assign phi
                phii[cellI] = recPhis[rCell];
            }
        }
    }
}

//---------------------------------------------------------------------------//
template <typename Type, typename volTypeField>
void ibInterpolation::interpolateToIntPoints
(
    volTypeField& phi,
    interpolation<Type>& interpPhi,
    List<List<intPoint>>& intInfos,
    List<List<Type>>& phiInIntPoints
)
{
    // prepare for sync
    List<DynamicList<Tuple2<label,label>>> labelsToSave(Pstream::nProcs());
    List<DynamicList<point>> iPointsToSync(Pstream::nProcs());
    List<DynamicList<label>> iCellsToSync(Pstream::nProcs());

    // loop over interpolation infos
    forAll(intInfos, iInfo)
    {
        // set size
        phiInIntPoints[iInfo].setSize(intInfos[iInfo].size());

        // loop over intPoints
        forAll(intInfos[iInfo], iPoint)
        {
            // Note (LK): think about the inclusion of surface point in interpolation points
            intPoint cPoint = intInfos[iInfo][iPoint];
            if (iPoint == 0) // skip the first intPoint (surface point)
            {
                phiInIntPoints[iInfo][iPoint] *= 0.0;
                continue;
            }

            // read data
            label iProc = cPoint.iProc_;

            if (Pstream::myProcNo() == iProc)
            {
                // interpolate and save
                Type phiP = interpPhi.interpolate(cPoint.iPoint_, cPoint.iCell_);
                phiInIntPoints[iInfo][iPoint] = phiP;
            }

            else
            {
                // save where to save
                Tuple2<label,label> labels;
                labels.first() = iInfo;
                labels.second() = iPoint;
                labelsToSave[iProc].append(labels);

                // save for sync
                iPointsToSync[iProc].append(cPoint.iPoint_);
                iCellsToSync[iProc].append(cPoint.iCell_);
            }
        }
    }

    // sync with other processors
    PstreamBuffers pBufsIPoints(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsICells(Pstream::commsTypes::nonBlocking);

    // send
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIPoints(proci, pBufsIPoints);
            UOPstream sendICells(proci, pBufsICells);
            sendIPoints << iPointsToSync[proci];
            sendICells << iCellsToSync[proci];
        }
    }

    pBufsIPoints.finishedSends();
    pBufsICells.finishedSends();

    // recieve
    List<DynamicList<Type>> phisToRetr(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIPoints(proci, pBufsIPoints);
            UIPstream recvICells(proci, pBufsICells);
            DynamicList<point> recIPoints (recvIPoints);
            DynamicList<label> recICells (recvICells);

            forAll(recIPoints, iPoint)
            {
                Type phiP = interpPhi.interpolate(recIPoints[iPoint], recICells[iPoint]);
                phisToRetr[proci].append(phiP);
            }
        }
    }

    // return
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendPhis(proci, pBufsIPoints);
            sendPhis << phisToRetr[proci];
        }
    }

    pBufsIPoints.finishedSends();

    List<DynamicList<Type>> phisCmpl(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvPhis(proci, pBufsIPoints);
            DynamicList<Type> recPhis (recvPhis);
            phisCmpl[proci] = recPhis;
        }
    }

    pBufsIPoints.clear();

    // save
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            forAll(labelsToSave[proci], iLabel)
            {
                // get the labels
                Tuple2<label,label> labels = labelsToSave[proci][iLabel];
                label iInfo = labels.first();
                label iPoint = labels.second();

                // save phi
                phiInIntPoints[iInfo][iPoint] = phisCmpl[proci][iLabel];
            }
        }
    }
}

// ************************************************************************* //
