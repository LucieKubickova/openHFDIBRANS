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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/
#include "lineIntInfo.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
lineIntInfo::lineIntInfo
(
    const fvMesh& mesh,
    ibMesh& ibMesh,
    List<label> ibCells,
    List<point> ibPoints,
    List<vector> ibNormals,
    bool correctIntPoints
)
:
mesh_(mesh),
ibMesh_(ibMesh),
ibCells_(ibCells),
ibPoints_(ibPoints),
ibNormals_(ibNormals),
correctIntPoints_(correctIntPoints)
{}
lineIntInfo::~lineIntInfo()
{}

//---------------------------------------------------------------------------//
void lineIntInfo::setIntpInfo
(
)
{
    // reset
    resetIntpInfo(ibCells_.size());

    // prepare interpolation points
    List<List<intPoint>>& intPoints = getIntPoints();

    // prepare lists
    List<DynamicList<point>> ibPointsToSolve(Pstream::nProcs()); // Note (LK): dunno if to include this as a part of the intPoint struct
    List<DynamicList<vector>> ibNormalsToSolve(Pstream::nProcs()); // Note (LK): dunno if to include this as a part of the intPoint struct
    List<DynamicList<intPoint>> intPointsToSolve(Pstream::nProcs());

    // save the first interpolation point 
    forAll (ibCells_, sCell)
    {
        // get surface cell label
        label cellI = ibCells_[sCell];
        //~ scalar intDist = Foam::pow(mesh_.V()[cellI],0.333); // Note (LK): int dist has to be bigger the 0.5*dX, but I am not sure about full dX, HFDIBDEM has full dX

        // save first interpolation point (surfPoint)
        intPoints[sCell].setSize(ORDER+1);
        intPoint cIntPoint
        (
            ibPoints_[sCell],
            cellI,
            Pstream::myProcNo(), // current processor
            Pstream::myProcNo(), // processor of origin (same for the first int point)
            sCell,               // label of the cell of origin
            false                // if its the last int point
        );
        intPoints[sCell][0] = cIntPoint;

        // save for looping lists
        ibPointsToSolve[Pstream::myProcNo()].append(ibPoints_[sCell]);
        ibNormalsToSolve[Pstream::myProcNo()].append(ibNormals_[sCell]);
        intPointsToSolve[Pstream::myProcNo()].append(cIntPoint);
    }

    // go by orders
    for(label i = 0; i < ORDER; ++i)
    {
        // lists to send
        List<DynamicList<point>> ibPointsToSend(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsToSend(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsToSend(Pstream::nProcs());

        // lists to continue
        List<DynamicList<point>> ibPointsToCont(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsToCont(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsToCont(Pstream::nProcs());

        // loop over processors
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            // loop over interpolation points
            forAll(intPointsToSolve[proci], iInfo)
            {
                // latest interpolation point
                intPoint cIntPoint = intPointsToSolve[proci][iInfo];
                point cPoint = cIntPoint.iPoint_;
                scalar intDist = Foam::pow(mesh_.V()[cIntPoint.iCell_],0.333); // Note (LK): int dist has to be bigger the 0.5*dX, but I am not sure about full dX, HFDIBDEM has full dX

                do {
                    cPoint += ibNormalsToSolve[proci][iInfo]*intDist;
                } while(ibMesh_.pointInCell(cPoint, cIntPoint.iCell_));

                // new interpolation points
                intPoint nIntPoint = findIntPoint(cIntPoint, cPoint);
                if (correctIntPoints_)
                {
                    correctIntPoint(ibPointsToSolve[proci][iInfo], nIntPoint);
                }

                // check for cells at domain boundary
                if (nIntPoint.iProc_ == -1)
                {
                    intPoints[cIntPoint.oLabel_][i].last_ = true;
                }

                // check if to send or keep
                else if (Pstream::myProcNo() != nIntPoint.iProc_)
                {
                    ibPointsToSend[nIntPoint.iProc_].append(ibPointsToSolve[proci][iInfo]);
                    ibNormalsToSend[nIntPoint.iProc_].append(ibNormalsToSolve[proci][iInfo]);
                    intPointsToSend[nIntPoint.iProc_].append(nIntPoint);
                }

                else
                {
                    ibPointsToCont[proci].append(ibPointsToSolve[proci][iInfo]);
                    ibNormalsToCont[proci].append(ibNormalsToSolve[proci][iInfo]);
                    intPointsToCont[proci].append(nIntPoint);
                }
            }
        }

        // sync with others
        List<DynamicList<point>> ibPointsRecv(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsRecv(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsRecv(Pstream::nProcs());
        sendAndRecvIntPoints(
                ibPointsToSend,
                ibNormalsToSend,
                intPointsToSend,
                ibPointsRecv,
                ibNormalsRecv,
                intPointsRecv);

        // clear lists
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            ibPointsToSend[proci].clear();
            ibNormalsToSend[proci].clear();
            intPointsToSend[proci].clear();
        }

        // finished solving of recieved points
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            forAll(intPointsRecv[proci], iInfo)
            {
                // get face label
                label faceI = intPointsRecv[proci][iInfo].iCell_;
                label cellI(0);
                
                // get cell label
                forAll(mesh_.boundaryMesh(), patchi)
                {
                    const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch
                            = refCast<const processorPolyPatch>(cPatch);

                        label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                            ? procPatch.neighbProcNo() : procPatch.myProcNo();

                        if (sProc == proci)
                        {
                            cellI = mesh_.faceOwner()[cPatch.start() + faceI];
                        }
                    }
                }

                // save
                intPointsRecv[proci][iInfo].iCell_ = cellI;

                // find interpolation point
                //~ point cPoint = mesh_.C()[cellI];
                intPoint cIntPoint
                (
                    //~ cPoint,
                    intPointsRecv[proci][iInfo].iPoint_,
                    cellI,
                    Pstream::myProcNo(),
                    intPointsRecv[proci][iInfo].oProc_,
                    intPointsRecv[proci][iInfo].oLabel_,
                    intPointsRecv[proci][iInfo].last_
                );

                //~ scalar intDist = Foam::pow(mesh_.V()[cellI],0.333); // Note (LK): int dist has to be bigger the 0.5*dX, but I am not sure about full dX, HFDIBDEM has full dX
                vector dir = cIntPoint.iPoint_ - ibPointsRecv[proci][iInfo];
                dir /= mag(dir);

                if (correctIntPoints_)
                {
                    correctIntPoint(ibPointsRecv[proci][iInfo], cIntPoint);
                }

                intPointsRecv[proci][iInfo] = cIntPoint;
            }
        }

        // return solved interpolation points to processor of origin
        List<DynamicList<intPoint>> intPointsSolved(Pstream::nProcs());
        returnSolvedIntPoints(intPointsToCont, intPointsRecv, intPointsSolved);
        
        // save solved points
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            forAll(intPointsSolved[proci], iInfo)
            {
                // get the original label
                label oLabel = intPointsSolved[proci][iInfo].oLabel_;

                // save int point
                intPoints[oLabel][i+1] = intPointsSolved[proci][iInfo];
            }
        }

        // clear lists and prepare for next order solution
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            ibPointsToSolve[proci].clear();
            ibNormalsToSolve[proci].clear();
            intPointsToSolve[proci].clear();

            forAll(intPointsToCont[proci], iInfo)
            {
                ibPointsToSolve[proci].append(ibPointsToCont[proci][iInfo]);
                ibNormalsToSolve[proci].append(ibNormalsToCont[proci][iInfo]);
                intPointsToSolve[proci].append(intPointsToCont[proci][iInfo]);
            }

            // clear lists
            ibPointsToCont[proci].clear();
            ibNormalsToCont[proci].clear();
            intPointsToCont[proci].clear();
        }

        // add recieved points
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            forAll(intPointsRecv[proci], iInfo)
            {
                ibPointsToSolve[proci].append(ibPointsRecv[proci][iInfo]);
                ibNormalsToSolve[proci].append(ibNormalsRecv[proci][iInfo]);
                intPointsToSolve[proci].append(intPointsRecv[proci][iInfo]);
            }

            // clear lists
            ibPointsRecv[proci].clear();
            ibNormalsRecv[proci].clear();
            intPointsRecv[proci].clear();
        }
    }
}

//---------------------------------------------------------------------------//
void lineIntInfo::correctIntPoint
(
    point ibPoint,
    intPoint& cPoint
)
{
    if(cPoint.iProc_ != Pstream::myProcNo())
    {
        return;
    }

    vector closestPoint = ibMesh_.getClosestPoint(ibPoint, cPoint);

    if(ibMesh_.pointInCell(closestPoint, cPoint.iCell_))
    {
        cPoint.iPoint_ = closestPoint;
    }
    else
    {
        const labelList& cellFaces(mesh_.cells()[cPoint.iCell_]);

        forAll (cellFaces, fi)
        {
            const face faceI = mesh_.faces()[cellFaces[fi]];
            vector dir = closestPoint - cPoint.iPoint_;

            pointHit pHit = faceI.ray(
                cPoint.iPoint_,
                dir,
                mesh_.points()
            );

            if(pHit.hit())
            {
                vector newP = 0.95*(pHit.hitPoint() - cPoint.iPoint_);
                newP += cPoint.iPoint_;

                if(ibMesh_.pointInCell(newP, cPoint.iCell_))
                {
                    cPoint.iPoint_ = newP;
                    break;
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
intPoint lineIntInfo::findIntPoint
(
    intPoint& fromP,
    point& endP
)
{
    if(fromP.iPoint_ == endP)
    {
        return intPoint();
    }

    intPoint retP
    (
        endP,
        fromP.iCell_,
        fromP.iProc_,
        fromP.oProc_,
        fromP.oLabel_,
        fromP.last_
    );

    if(fromP.iProc_ == Pstream::myProcNo())
    {
        label faceInDir = -1;
        while(!ibMesh_.pointInCell(retP.iPoint_, retP.iCell_))
        {
            //~ faceInDir = ibMesh_.getFaceInDir(retP, faceInDir);
            vector dir = retP.iPoint_ - mesh_.C()[retP.iCell_];
            faceInDir = ibMesh_.getFaceInDir(retP.iCell_, dir);
            if (!mesh_.isInternalFace(faceInDir))
            {
                label facePatchId(mesh_.boundaryMesh().whichPatch(faceInDir));
                const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];

                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch
                        = refCast<const processorPolyPatch>(cPatch);
                    label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                        ? procPatch.neighbProcNo() : procPatch.myProcNo();

                    retP.iCell_ = cPatch.whichFace(faceInDir);
                    retP.iProc_ = sProc;

                    return retP;
                }
                else
                {
                    retP.iProc_ = -1;
                    return retP;
                }
            }

            label owner(mesh_.owner()[faceInDir]);
            label neighbour(mesh_.neighbour()[faceInDir]);
            retP.iCell_ = (retP.iCell_ == neighbour) ? owner : neighbour;
        }

        return retP;
    }

    return retP;
}

//---------------------------------------------------------------------------//
void lineIntInfo::sendAndRecvIntPoints
(
    List<DynamicList<point>>& ibPointsToSend,
    List<DynamicList<vector>>& ibNormalsToSend,
    List<DynamicList<intPoint>>& intPointsToSend,
    List<DynamicList<point>>& ibPointsRecv,
    List<DynamicList<vector>>& ibNormalsRecv,
    List<DynamicList<intPoint>>& intPointsRecv
)
{
    PstreamBuffers pBufsIbPoints(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsIbNormals(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsIntPoints(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream sendIbPoints(proci, pBufsIbPoints);
            UOPstream sendIbNormals(proci, pBufsIbNormals);
            UOPstream sendIntPoints(proci, pBufsIntPoints);
            sendIbPoints << ibPointsToSend[proci];
            sendIbNormals << ibNormalsToSend[proci];
            sendIntPoints << intPointsToSend[proci];
        }
    }

    pBufsIbPoints.finishedSends();
    pBufsIbNormals.finishedSends();
    pBufsIntPoints.finishedSends();

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIbPoints(proci, pBufsIbPoints);
            UIPstream recvIbNormals(proci, pBufsIbNormals);
            UIPstream recvIntPoints(proci, pBufsIntPoints);
            DynamicList<point> recIbPoints (recvIbPoints);
            DynamicList<point> recIbNormals (recvIbNormals);
            DynamicList<intPoint> recIntPoints (recvIntPoints);
            ibPointsRecv[proci] = recIbPoints;
            ibNormalsRecv[proci] = recIbNormals;
            intPointsRecv[proci] = recIntPoints;
        }
    }

    pBufsIbPoints.clear();
    pBufsIbNormals.clear();
    pBufsIntPoints.clear();
}

//---------------------------------------------------------------------------//
void lineIntInfo::returnSolvedIntPoints
(
    List<DynamicList<intPoint>>& intPointsToCont,
    List<DynamicList<intPoint>>& intPointsRecv,
    List<DynamicList<intPoint>>& intPointsSolved
)
{
    // prepare list to return
    List<DynamicList<intPoint>> intPointsToRetr(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        forAll(intPointsToCont[proci], iInfo)
        {
            // get the label of the processor of origin
            label oProc = intPointsToCont[proci][iInfo].oProc_;

            // save to solved 
            if (Pstream::myProcNo() == oProc)
            {
                intPointsSolved[oProc].append(intPointsToCont[proci][iInfo]);
            }

            // add to return
            else
            {
                intPointsToRetr[oProc].append(intPointsToCont[proci][iInfo]);
            }
        }

        forAll(intPointsRecv[proci], iInfo)
        {
            // get the label of the processor of origin
            label oProc = intPointsRecv[proci][iInfo].oProc_;

            // save to solved 
            if (Pstream::myProcNo() == oProc)
            {
                intPointsSolved[oProc].append(intPointsRecv[proci][iInfo]);
            }

            // add to return
            else
            {
                intPointsToRetr[oProc].append(intPointsRecv[proci][iInfo]);
            }
        }
    }

    // send
    PstreamBuffers pBufsIntPoints(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream sendIntPoints(proci, pBufsIntPoints);
            sendIntPoints << intPointsToRetr[proci];
        }
    }

    pBufsIntPoints.finishedSends();

    // recieve
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIntPoints(proci, pBufsIntPoints);
            DynamicList<intPoint> recIntPoints (recvIntPoints);
            intPointsSolved[proci] = recIntPoints;
        }
    }

    pBufsIntPoints.clear();
}

//---------------------------------------------------------------------------//
