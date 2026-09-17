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
    openHFDIBRANS is licensed under the GNU LESSER GENERAL PUBLIC LICENSE
    (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the
    terms and conditions of version 3 of the GNU General Public License,
    supplemented by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIBRANS. If not, see
    <http://www.gnu.org/licenses/lgpl.html>.

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),

\*---------------------------------------------------------------------------*/

#include "ibCutCell.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ibCutCell::~ibCutCell()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void ibCutCell::makeFaceCentresAndAreas()
{
    if (type_ != 0)
	{
		return;
	};

    faceAreas_.setSize(faces_.size());
    faceCentres_.setSize(faces_.size());

    forAll (faces_, facei)
    {
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
                const vector thisPoint = f[pi];
                const vector nextPoint = f[(pi + 1)%nPoints];

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

// ------------------------------------------------------------------------- //

void ibCutCell::syncNormals()
{
    if (type_ != 0)
	{
		return;
	};

    const DynamicList<DynamicList<point>>& cFaces = faces_;

    vector C = C_;

    forAll (cFaces, facei)
    {
        const point& fc = faceCentres_[facei];
        vector d = fc - C;

        if
		(
			(d/(mag(d) + VSMALL)
		  & faceAreas_[facei]/(mag(faceAreas_[facei]) + VSMALL)) < -SMALL
		)
		{
			faceAreas_[facei] *= -1.0;
		};
    }

}

// ------------------------------------------------------------------------- //

void ibCutCell::makeCellCentreAndVol()
{
    if (type_ != 0)
	{
		return;
	};

    const DynamicList<DynamicList<point>>& cFaces = faces_;

    // Restart possibly wrong data
    C_ *= 0.0;
    V_ *= 0.0;

    // Estimate the cell centre and bounding box using the face centres
    vector cEst = Zero;

    forAll (cFaces, facei)
    {
        const point& fc = faceCentres_[facei];
        cEst += fc;
    }
    cEst /= cFaces.size();


    // Sum up the face-pyramid contributions
    forAll (cFaces, facei)
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar ibCutCell::yOrtho()
{
    const labelList& inLabels = inCell_.labels(mesh_.faces());
    const List<point>& inPts = inCell_.points(mesh_.faces(),mesh_.points());
    const edgeList& inEdges = inCell_.edges(mesh_.faces());

    label nPoints = inLabels.size();

    boolList isIn(nPoints, true);
    label nVertsOut = 0;
    forAll(isIn, bi)
    {
        if (!(((inPts[bi]-surfPoint_) & surfNorm_) < -SMALL))
        {
            isIn[bi] = false;
            nVertsOut += 1;
        }

    }

    if (nVertsOut == 0)
    {
        type_ = -1;

        V_ = VIn_;
        C_ = CIn_;
    }
    if (nVertsOut == nPoints)
    {
        type_ = 1;
        V_ = VIn_;
        C_ = CIn_;
    }

    bool isCut = (nVertsOut > 0) && (nVertsOut < nPoints);

    // Note (MI): if all is inside, I should just reconstruct the
    //            current cell
    if (isCut)
    {
        type_ = 0;

        // Construct all new faces but the purely cutFace
        forAll (inCell_, facei)
        {
            const face& f = mesh_.faces()[inCell_[facei]];

            const List<point>& fPoints = f.points(mesh_.points());
            label nFPoints = fPoints.size();
            label nVertsOut = 0;

            DynamicList<point> fV;

            forAll(fPoints, fpi)
            {
                point sP = fPoints[fpi];
                point eP = fPoints[(fpi + 1) % nFPoints];

                scalar sdS = (sP - surfPoint_) & surfNorm_;
                scalar sdE = (eP - surfPoint_) & surfNorm_;

                if (sign(sdS/mag(sP - surfPoint_)) > 0.0)
                {
                    nVertsOut += 1;
                    fV.append(sP);
                }

                if
				(
					sign(sdS/mag(sP - surfPoint_)*sdE/mag(eP - surfPoint_))
				 <= 0.0
				)
                {
                    vector eVec = (eP - sP);
                    eVec /= mag(eVec);

                    point cP =
					(
						sP
					  + ((surfPoint_ - sP) & surfNorm_)
					   /(eVec & surfNorm_)
					   *eVec
					);

                    // cutCell will have cP instead of sP
                    if (mag(cP - sP) < 2.0*mag(eP - sP))
                    {
                        fV.append(cP);
                    }
                    else
                    {
                        fV.append(sP);
                    }
                }
            }

			// Note (VV): Ask Lucka about this *unique* if statement
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

        // Construct the cutFace
        DynamicList<point> fV;
        label cutFaces = 0;
        label gFacei = -1;

        List<label> unCheckedFaces = inCell_;
        DynamicList<label> checkedEdges;

        forAll(inCell_, facei)
        {
            if (cutFaces == 0)
			{
				gFacei = inCell_[facei];
			}

            forAll(unCheckedFaces, ufi)
            {
                if (unCheckedFaces[ufi] == gFacei)
                {
                    unCheckedFaces[ufi] = -1;
                    break;
                }
            }

            const labelList& fEdges = mesh_.faceEdges()[gFacei];

            label nCutEdges = 0;
            label cutEdge = -1;

            forAll(fEdges,edgei)
            {
                const edge& e = mesh_.edges()[fEdges[edgei]];

                point sP = mesh_.points()[e.start()];
                point eP = mesh_.points()[e.end()];

                scalar sdS((sP - surfPoint_) & surfNorm_);
                scalar sdE((eP - surfPoint_) & surfNorm_);

                if (sign(sdS*sdE) < 0.0 or sdS == 0.0)
                {
                    bool toInclude = true;

                    forAll(checkedEdges, cei)
                    {
                        if (checkedEdges[cei] == fEdges[edgei])
                        {
                            toInclude = false;
                            break;
                        }
                    }

                    vector eVec = eP - sP;
                    eVec /= mag(eVec);

                    point cP = vector::zero;
                    if (sdS == 0.0)
                    {
                        cP = sP;
                    }
                    else
                    {
                        cP =
						(
							sP
						  + ((surfPoint_ - sP) & surfNorm_)
						   /(eVec & surfNorm_)
						   *eVec
						);
                    }

                    // cutCell will have cP instead of sP
                    if (toInclude)
                    {
                        cutEdge = fEdges[edgei];
                        // Note (MI): this is important if I first find
                        //            an unchecked cutEdge and afterwards
                        //            the checked one
                        // Note (MI): each cutFace has 2 cutEdges

                        if (mag(cP - sP) < 2.0*mag(eP - sP))
                        {
                            fV.append(cP);
                        }
                        else
                        {
                            fV.append(sP);
                        }
                        checkedEdges.append(cutEdge);
                    }
                    nCutEdges++;
                }

                if (nCutEdges == 2 && cutEdge > 0)
                {
                    // Get the last found cutEdge faces
                    const labelList& edgeFaces = mesh_.edgeFaces()[cutEdge];

                    // Identify face neighbor in the current cell
                    forAll(edgeFaces,eFacei)
                    {
                        if (edgeFaces[eFacei] != gFacei)
                        {
                            bool toBreak = false;
                            forAll(unCheckedFaces,cFacei)
                            {
                                if
								(
									edgeFaces[eFacei] == unCheckedFaces[cFacei]
								)
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

            if (cutFaces == nCutFaces_)
			{
				break;
			}
        }
        v_.append(fV);
    }
    faces_ = v_;
	// Note (MI): in my v_, I have basically a face structure, but
    //            my face structure does not contain labels to individual
    //            vertices. Instead, it contains the complete data

    // Note (LK): do not bother with uncut cells
    if (isCut)
    {
        makeFaceCentresAndAreas();

        makeCellCentreAndVol();

        syncNormals();
    }

    return ((C_-surfPoint_) & surfNorm_);
}

// ------------------------------------------------------------------------- //

scalar ibCutCell::yOrthoEst()
{
    scalar dotMax = 0.0;
    scalar yO = 0.0;

    forAll (inCell_, facei)
    {
        const face& f = mesh_.faces()[inCell_[facei]];

        point fC = f.centre(mesh_.points());
        vector fN = f.normal(mesh_.points());

		// Flip normals if required
        if ((fN & (fC - CIn_)) < -SMALL)
		{
			fN *= -1.0;
		};

        scalar d = (surfPoint_ - fC) & fN;

        scalar dotC = fN & surfNorm_;
        if (dotC > dotMax)
        {
            dotMax = dotC;
            yO     = d;
        }
    }

    type_ = 0;

	// If this fails, use the costly evaluation
    if (yO*(-0.5) < SMALL)
    {
        Info<< "Failsafe to yOrtho" << endl;
        yO = yOrtho()*(-2.0);
    };

    return max(yO*(-0.5), SMALL);
}


// ************************************************************************* //
