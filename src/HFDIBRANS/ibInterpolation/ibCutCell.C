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
    cutCell class to be used for hfdibRANS turbulence models

Contributors
    Martin Isoz (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/
#include "ibCutCell.H"

#include "ListOps.H"

using namespace Foam;

ibCutCell::ibCutCell
(
    const fvMesh& mesh,
    const vector& surfNorm,
    const point& surfPoint,
    const cell& inCell
)
:
mesh_(mesh),
surfNorm_(surfNorm),
surfPoint_(surfPoint),
inCell_(inCell),
type_(-2),
nCutFaces_(0)
{
    CIn_ = inCell_.centre(mesh_.points(),mesh_.faces());
    VIn_ = inCell_.mag(mesh_.points(),mesh_.faces());
}

//---------------------------------------------------------------------------//
ibCutCell::~ibCutCell()
{
}

//---------------------------------------------------------------------------//
void ibCutCell::makeFaceCentresAndAreas()
{
    if (type_ != 0) {return;};
    
    faceAreas_.setSize(faces_.size());
    faceCentres_.setSize(faces_.size());
    
    forAll (faces_,facei)
    {
        // const face& f = faces_[facei];
        const DynamicList<point>& f = faces_[facei];
        const label nPoints = f.size();
         
        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            faceCentres_[facei] = (1.0/3.0)*(f[0] + f[1] + f[2]);
            faceAreas_[facei] = 0.5*((f[1] - f[0])^(f[2] - f[0]));
        }
        else
        {
            vector sumN = Zero;
            scalar sumA = Zero;
            vector sumAc = Zero;
 
            vector fCentre = f[0];
            for (label pi = 1; pi < nPoints; ++pi)
            {
                fCentre += vector(f[pi]);
            }
            fCentre /= nPoints;
 
            for (label pi = 0; pi < nPoints; ++pi)
            {
                const vector thisPoint(f[pi]);
                const vector nextPoint(f[(pi+1)%nPoints]);
 
                vector c = thisPoint + nextPoint + fCentre;
                vector n = (nextPoint - thisPoint)^(fCentre - thisPoint);
                scalar a = mag(n);
 
                sumN += n;
                sumA += a;
                sumAc += a*c;
            }
 
            if (sumA < ROOTVSMALL)
            {
                faceCentres_[facei] = fCentre;
                faceAreas_[facei] = Zero;
            }
            else
            {
                faceCentres_[facei] = (1.0/3.0)*sumAc/sumA;
                faceAreas_[facei] = 0.5*sumN;
            }
        }
    }
}

void ibCutCell::syncNormals()
{
    if (type_ != 0) {return;};
    
    const DynamicList<DynamicList<point>>& cFaces = faces_;
    
    vector C(C_);
    
    forAll (cFaces,facei)
    {
        const point& fc = faceCentres_[facei];
        vector d(fc - C);
        
        if ( (d/(mag(d)+VSMALL) & faceAreas_[facei]/(mag(faceAreas_[facei])+VSMALL)) < -SMALL) {faceAreas_[facei]*=-1.0;};
    }
    
}

void ibCutCell::makeCellCentreAndVol()
{
    if (type_ != 0) {return;};
    
    const DynamicList<DynamicList<point>>& cFaces = faces_;
    
    // restart possibly wrong data
    C_ *= 0.0;
    V_ *= 0.0;
 
    // Estimate the cell centre and bounding box using the face centres
    vector cEst(Zero);

    forAll (cFaces,facei)
    {
        const point& fc = faceCentres_[facei];
        cEst += fc;
    }
    cEst /= cFaces.size();


    // Sum up the face-pyramid contributions
    forAll (cFaces,facei)
    {
        
        // Calculate 3* the face-pyramid volume
        scalar pyr3Vol = mag(faceAreas_[facei] & (faceCentres_[facei] - cEst));

        // Accumulate face-pyramid volume
        V_ += pyr3Vol;

        // Calculate the face-pyramid centre
        const vector pCtr = (3.0/4.0)*faceCentres_[facei] + (1.0/4.0)*cEst;
        
        // Accumulate volume-weighted face-pyramid centre
        C_ += pyr3Vol*pCtr;
    }

    // Average the accumulated quantities

    if (mag(V_) > VSMALL)
    {
        point cc = C_ / V_;
        
        C_ = cc;
    }
    else
    {
        C_ = cEst;
    }

    V_ *= (1.0/3.0);
}

scalar ibCutCell::yOrtho()
{    
    const labelList& inLabels(inCell_.labels(mesh_.faces()));
    const List<point>& inPts(inCell_.points(mesh_.faces(),mesh_.points()));
    const edgeList& inEdges(inCell_.edges(mesh_.faces()));
    
    label nPoints(inLabels.size());
    
    boolList isIn(nPoints,true);
    label nVertsOut(0);
    forAll (isIn,bi)
    {
        if (((inPts[bi]-surfPoint_) & surfNorm_) > SMALL)
        {
            isIn[bi] = false;
            nVertsOut+= 1;
        }
        
    }
    
    if (nVertsOut == 0) 
    {
        //~ Info << "cell is in" << endl;
        type_ = -1;
        
        V_ = VIn_;
        C_ = CIn_;
    };
    if (nVertsOut == nPoints)
    {
        //~ Info << "cell is out" << endl;
        type_ = 1;
        V_ = VIn_;
        C_ = CIn_;
    };
    
    bool isCut((nVertsOut > 0) and (nVertsOut < nPoints));
    // Note (MI): if all is inside, I should just reconstruct the
    //            current cell
    if (isCut)
    {            
        type_ = 0;
        
        // construct all new faces but the purely cutFace
        forAll (inCell_,facei)
        {
            const face& f = mesh_.faces()[inCell_[facei]];
            
            const List<point>& fPoints(f.points(mesh_.points()));
            label nFPoints(fPoints.size());
            label nVertsOut(0);
                
            DynamicList<point> fV;
            
            forAll(fPoints,fpi)
            {
                point sP(fPoints[fpi]);
                point eP(fPoints[(fpi+1) % nFPoints]);
                
                scalar sdS((sP-surfPoint_) & surfNorm_);
                scalar sdE((eP-surfPoint_) & surfNorm_);

                if (sign(sdS/mag(sP-surfPoint_)) > SMALL)
                {
                    nVertsOut+= 1;
                    fV.append(sP);
                }
                else if (sign(sdS/mag(sP-surfPoint_) * sdE/mag(eP-surfPoint_)) < 0.0)
                {
                    
                    vector eVec((eP - sP));
                    eVec /= mag(eVec);
                    
                    point cP = sP + ((surfPoint_ - sP) & surfNorm_) / (eVec & surfNorm_)*eVec;
                    
                    // cutCell will have cP instead of sP
                    
                    if (mag(cP-sP) < 2.0*mag(eP-sP))
                    {
                        fV.append(cP);
                    }
                    else
                    {
                        fV.append(sP);
                    }
                }
            }
            
            if (nVertsOut == 0)
            {
                //~ Info << "fully inner face" << endl;
            }            
            else
            {
                v_.append(fV);
                if (nVertsOut == nFPoints)
                {
                    //~ Info << "fully outer face" << endl;
                }
                else
                {
                    nCutFaces_++;
                }
            }
        }
                
        // construct the cutFace
        DynamicList<point> fV;                                          //list to contain cutPoints forming cutFace
        label cutFaces(0);                                              //cutFaces counter
        label gFacei(-1);                                               //examined cutFace index
        
        List<label> unCheckedFaces(inCell_);                            //initialized yet unchecked cell faces
        DynamicList<label> checkedEdges;                                //initialize already checked edges
        
        forAll (inCell_,facei)                                          //for all input cell faces
        {
            if (cutFaces == 0) {gFacei = inCell_[facei];}               //initialization via input face index
            
            forAll (unCheckedFaces, ufi)                                //replace face index in unCheckedFaces by -1
            {
                if (unCheckedFaces[ufi] == gFacei)
                {
                    unCheckedFaces[ufi] = -1;
                    break;
                }
            }
                        
            const labelList& fEdges = mesh_.faceEdges()[gFacei];        //get labels of face edges
            
            label nCutEdges(0);                                         //initialize cutEdges counter
            label cutEdge(-1);                                          //see if the current edge is a new cutEdge
            
            forAll  (fEdges,edgei)                                      //examine the face edges
            {
                const edge& e = mesh_.edges()[fEdges[edgei]];
                                
                point sP(mesh_.points()[e.start()]);
                point eP(mesh_.points()[e.end()]);
                
                scalar sdS((sP-surfPoint_) & surfNorm_);
                scalar sdE((eP-surfPoint_) & surfNorm_);
                
                if (sign(sdS * sdE) < 0.0)                              //found a cut edge
                {                    
                    bool toInclude(true);                               //I want to include this
                                        
                    forAll (checkedEdges, cei)                          //if the edge is already included, do not include
                    {
                        if (checkedEdges[cei] == fEdges[edgei])
                        {
                            toInclude = false;
                            break;
                        }
                    }
                    
                    vector eVec((eP - sP));
                    eVec /= mag(eVec);
                    
                    point cP = sP + ((surfPoint_ - sP) & surfNorm_) / (eVec & surfNorm_)*eVec;
                    
                    // cutCell will have cP instead of sP
                    
                    if (toInclude)                                      //yet not included cutEdge
                    {
                        cutEdge = fEdges[edgei];                        //save the edge index
                        // Note (MI): this is important if I first find
                        //            an unchecked cutEdge and afterwards
                        //            the checked one
                        // Note (MI): yeach cutFace has 2 cutEdges
                        
                        if (mag(cP-sP) < 2.0*mag(eP-sP))                //append only if the cutPoint is "close" to the current edge
                        {
                            fV.append(cP);
                        }
                        else
                        {
                            fV.append(sP);
                        }
                        checkedEdges.append(cutEdge);                   //include in the already checked edges
                    }
                    nCutEdges++;
                }
                
                if (nCutEdges == 2 and cutEdge > 0)                     //found two cut edges on the current face, finish with the current face
                {
                    // get the last found cutEdge faces
                    const labelList& edgeFaces = mesh_.edgeFaces()[cutEdge];
                    
                    // identify face neighbor in the current cell
                    forAll (edgeFaces,eFacei)
                    {
                        if (edgeFaces[eFacei] != gFacei)
                        {
                            bool toBreak(false);
                            forAll(unCheckedFaces,cFacei)
                            {
                                if (edgeFaces[eFacei] == unCheckedFaces[cFacei])
                                {
                                    gFacei = unCheckedFaces[cFacei];
                                    toBreak = true;
                                    break;
                                }
                            }
                            if (toBreak) {break;}
                        }
                    }
                    
                    cutFaces++;
                    break;
                }
                
            }
            
            if (cutFaces == nCutFaces_) {break;}                        //I examined all the cutFaces
            
        }
        
        v_.append(fV);
    }
    faces_ = v_;                                                        //ugly, ugly, ugly
    // Note (MI): in my v_, I have basically a face structure, but
    //            my face structure does not contain labels to individual
    //            vertices. Instead, it contains the complete data
    
    makeFaceCentresAndAreas();
        
    makeCellCentreAndVol();
    
    syncNormals();
        
    return ((C_-surfPoint_) & surfNorm_);
}

scalar ibCutCell::yOrthoEst()
{
    scalar dotMax(0.0);
    scalar yO(0.0);
    
    forAll (inCell_,facei)
    {
        const face& f = mesh_.faces()[inCell_[facei]];
        
        point   fC(f.centre(mesh_.points()));
        vector  fN(f.normal(mesh_.points()));
        
        if ((fN & (fC - CIn_)) < -SMALL) {fN *= -1.0;};                 //flip normals if required
        
        scalar d = (surfPoint_ - fC) & fN;
        
        scalar dotC(fN & surfNorm_);
        if (dotC > dotMax)
        {
            dotMax = dotC;
            yO     = d;
        }
    }
    
    type_ = 0;                                                          //automatically assumed cut
    
    if (yO*(-0.5) < SMALL)                                              //if it fails, use the costly value
    {
        Info << "failsafe to yOrtho" << endl;
        yO = yOrtho()*(-2.0);
    };
    
    return max(yO*(-0.5),SMALL);
}
