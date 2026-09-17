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

#include "ibDirichletBCs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ibDirichletBCs::ibDirichletBCs
(
    const fvMesh& mesh,
    ibMesh& ibMesh,
    const volScalarField& body,
    List<DynamicList<boundaryCell>>& boundaryCells,
    List<DynamicList<surfaceCell>>& surfaceCells,
    labelField& isBoundaryCell
)
:
	mesh_(mesh),
	ibMesh_(ibMesh),
	body_(body),
	boundaryCells_(boundaryCells),
	surfaceCells_(surfaceCells),
	isBoundaryCell_(isBoundaryCell),
	turbulenceProperties_
	(
		IOobject
		(
			"turbulenceProperties",
			"constant",
			mesh_,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	),
	HFDIBDEMDict_
	(
		IOobject
		(
			"HFDIBDEMDict",
			"constant",
			mesh_,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	),
	fvSchemes_
	(
		IOobject
		(
			"fvSchemes",
			"system",
			mesh_,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	),
	yPlusi_
	(
		IOobject
		(
			"yPlusi",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedScalar("zero", dimless, -1.0)
	),
	uTaui_
	(
		IOobject
		(
			"uTaui",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedScalar("zero", dimless, 0.0)
	),
	nuti_
	(
		IOobject
		(
			"nuti",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedScalar("zero", dimless, 0.0)
	),
	kappa_(0.41),
	E_(9.8),
	Cmu_(0.09),
	Ceps2_(1.9),
	beta1_(0.075)
{
    // Initiate lists
    nutAtIB_.setSize(Pstream::nProcs());
    kAtIB_.setSize(Pstream::nProcs());
    uTauAtIB_.setSize(Pstream::nProcs());

    // Read turbulence properties
    turbulenceProperties_.lookup("simulationType") >> simulationType_;

    // Read HFDIBDEMDict
    thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
    useYEff_ = HFDIBDEMDict_.lookupOrDefault<bool>("useEffectiveDist", true);
    uTauType_ = HFDIBDEMDict_.lookupOrDefault<word>("uTauType", "freeStreamCell");
    uTauCoeff_ = HFDIBDEMDict_.lookupOrDefault<scalar>("uTauCoeff", 1.0);

    // Read boundary condition for velocity
    HFDIBBCsDict_ = HFDIBDEMDict_.subDict("wallFunctions");
    UBC_ = HFDIBBCsDict_.lookupOrDefault<word>("U", "noSlip");

    // Read simulation type
    if (simulationType_ != "laminar")
    {
        HFDIBBCsDict_.lookup("nut") >> nutWF_;
        HFDIBBCsDict_.lookup("k") >> kWF_;
        HFDIBBCsDict_.lookup("omega") >> omegaWF_;
        HFDIBBCsDict_.lookup("epsilon") >> epsilonWF_;
    }

    // Compute turbulence parameters
    Cmu75_ = Foam::pow(Cmu_, 0.75);
    Cmu25_ = pow025(Cmu_);
    Cmu5_ = Foam::sqrt(Cmu_);
    calcYPlusLam();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ibDirichletBCs::~ibDirichletBCs()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void ibDirichletBCs::calcYPlusLam()
{
    yPlusLam_ = 11.0;

    for (int i = 0; i < 10; i++)
    {
        yPlusLam_ = Foam::log(max(E_*yPlusLam_, 1))/kappa_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ibDirichletBCs::setSizeToLists()
{
    // Set size
    nutAtIB_[Pstream::myProcNo()].setSize
	(
		boundaryCells_[Pstream::myProcNo()].size()
	);
    kAtIB_[Pstream::myProcNo()].setSize
	(
		boundaryCells_[Pstream::myProcNo()].size()
	);
    uTauAtIB_[Pstream::myProcNo()].setSize
	(
		boundaryCells_[Pstream::myProcNo()].size()
	);

    // Reset
    nutAtIB_[Pstream::myProcNo()] = 0.0;
    kAtIB_[Pstream::myProcNo()] = 0.0;
    uTauAtIB_[Pstream::myProcNo()] = 0.0;
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::UAtIB
(
    List<vector>& UIB,
    volVectorField& U
)
{
    if (simulationType_ == "laminar" or UBC_ == "noSlip")
    {
        forAll(UIB, uCell)
        {
            // Assign zero
            UIB[uCell] = ibZero(UIB[uCell]);
        }
    }
    else if (UBC_ == "partialSlip")
    {
        // Read partial clip coefficient
        scalar alpha = HFDIBBCsDict_.lookupOrDefault<scalar>("UCoeff", 0.0);

        // Loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Reset field
            UIB[bCell] = vector::zero;

            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // Get surface normal
            vector normal = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;

            // Calculate value at boundary
            UIB[bCell] += (1 - alpha)*transform(I - sqr(normal), U[cellI]);
        }
    }
    else
    {
        FatalError
			<< UBC_ << " condition for U in " << simulationType_
			<< " not implemented at the IB" << exit(FatalError);
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::TAtIB
(
    List<scalar>& TIB,
    scalar TIn
)
{
    forAll(TIB, bCell)
    {
        // Assign value
        TIB[bCell] = TIn;
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::updateUTauAtIB
(
    volScalarField& k
)
{
    // Prepare interpolation scheme
    dictionary HFDIBInnerSchemes =
		fvSchemes_.subDict("HFDIBSchemes").subDict("innerSchemes");
    autoPtr<interpolation<scalar>> interpK =
		interpolation<scalar>::New(HFDIBInnerSchemes, k);

    // Prepare synchronization
    List<DynamicList<label>> fCellsToSync(Pstream::nProcs());
    List<DynamicList<point>> fPointsToSync(Pstream::nProcs());
    List<DynamicList<label>> bLabelsToRecv(Pstream::nProcs());

    // If uTau from boundary cell
    if (uTauType_ == "boundaryCell")
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Reset field
            uTauAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // Compute uTau
            uTauAtIB_[Pstream::myProcNo()][bCell] =
				Cmu25_*Foam::sqrt(k[cellI]);
        }
    }
    else if
	(
		uTauType_ == "effectiveDistance"
	 || uTauType_ == "cellDistance"
	 || uTauType_ == "coeffDistance"
	)
    {
        // Loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Reset field
            uTauAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
            label fCell1 = boundaryCells_[Pstream::myProcNo()][bCell].fCell1_;
            label fCell2 = boundaryCells_[Pstream::myProcNo()][bCell].fCell2_;
            label fProc1 = boundaryCells_[Pstream::myProcNo()][bCell].fProc1_;
            label fProc2 = boundaryCells_[Pstream::myProcNo()][bCell].fProc2_;

            // Get surface point and normal
            point sPoint = boundaryCells_[Pstream::myProcNo()][bCell].sPoint_;
            vector sNorm = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;

            // Get distance
            scalar dist = 0.0;
            if (uTauType_ == "effectiveDistance")
            {
                dist = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else if (uTauType_ == "cellDistance")
            {
                dist = ibMesh_.getCellSize(cellI, sNorm);
            }
            else if (uTauType_ == "coeffDistance")
            {
                dist = 1.0; // Note (LK): value set by uTauCoeff
            }

            // fictional point in effective distance
            point distPoint = sPoint + uTauCoeff_*dist*sNorm;

            // check in which cell the point is
            label kCell;
            label kProc;
            if (ibMesh_.pointInCell(distPoint, cellI))
            {
                kCell = cellI;
                kProc = Pstream::myProcNo();
            }
			//~ else if (ibMesh_.pointInCell(yEffPoint, fCell1)) // Note (LK): the distPoint should not be farther than one cell away
            else
            {
                kCell = fCell1;
                kProc = fProc1;
            }
            //~ else
            //~ {
                //~ kCell = fCell2;
                //~ kProc = fProc2;
            //~ }

            if (kProc == Pstream::myProcNo())
            {
                // Interpolate k
                scalar kPoint = interpK->interpolate(distPoint, kCell);

                // Compute friction velocity
                uTauAtIB_[Pstream::myProcNo()][bCell] =
					Cmu25_*Foam::sqrt(kPoint);
            }
            else
            {
                fCellsToSync[kProc].append(kCell);
                fPointsToSync[kProc].append(distPoint);
                bLabelsToRecv[kProc].append(bCell);
            }
        }
    }
    else if (uTauType_ == "freeStreamCell" || uTauType_ == "interpPoint")
    {
        // Loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Reset field
            uTauAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // Prepare variables
            label fCell = 0;
            label fProc = -1;
            point fPoint = vector::zero;

            // Use the first interpolation point
            fCell = boundaryCells_[Pstream::myProcNo()][bCell].fCell1_;
            fProc = boundaryCells_[Pstream::myProcNo()][bCell].fProc1_;
            fPoint = boundaryCells_[Pstream::myProcNo()][bCell].fPoint1_;

            // Check body -- Note (LK): not really working for pipe or anything, but idea to smooth out
            //~ if (body_[cellI] < thrSurf_)
            //~ {
                //~ fCell = boundaryCells_[Pstream::myProcNo()][bCell].fCell1_;
                //~ fProc = boundaryCells_[Pstream::myProcNo()][bCell].fProc1_;
                //~ fPoint = boundaryCells_[Pstream::myProcNo()][bCell].fPoint1_;
            //~ }
            //~ else
            //~ {
                //~ fCell = boundaryCells_[Pstream::myProcNo()][bCell].fCell2_;
                //~ fProc = boundaryCells_[Pstream::myProcNo()][bCell].fProc2_;
                //~ fPoint = boundaryCells_[Pstream::myProcNo()][bCell].fPoint2_;
            //~ }

            // Compute uTau based on values from the free stream
            if (Pstream::myProcNo() == fProc)
            {
                if (uTauType_ == "interpPoint")
                {
                    // Interpolate k
                    scalar kPoint = interpK->interpolate(fPoint, fCell);

                    // Compute friction velocity
                    uTauAtIB_[Pstream::myProcNo()][bCell] =
						Cmu25_*Foam::sqrt(kPoint);
                }
                else
                {
                    // Compute friction velocity
                    uTauAtIB_[Pstream::myProcNo()][bCell] =
						Cmu25_*Foam::sqrt(k[fCell]);
                }
            }
            else if (fProc == -1)
            {
                // Get uTau from the boundary cell itself
                uTauAtIB_[Pstream::myProcNo()][bCell] =
					Cmu25_*Foam::sqrt(k[cellI]);
            }
            else
            {
                fCellsToSync[fProc].append(fCell);
                fPointsToSync[fProc].append(fPoint);
                bLabelsToRecv[fProc].append(bCell);
            }
        }
    }
    else
    {
        FatalError
			<< "uTau calculation type " << uTauType_
			<< " not implemented" << exit(FatalError);
    }

    // Sync with other processors
    PstreamBuffers pBufsFCells(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsFPoints(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendFCells(proci, pBufsFCells);
            UOPstream sendFPoints(proci, pBufsFPoints);
            sendFCells << fCellsToSync[proci];
            sendFPoints << fPointsToSync[proci];
        }
    }

    pBufsFCells.finishedSends();
    pBufsFPoints.finishedSends();

    // Recieve
    List<DynamicList<scalar>> uTausToRetr(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvFCells(proci, pBufsFCells);
            UIPstream recvFPoints(proci, pBufsFPoints);
            DynamicList<label> recFCells (recvFCells);
            DynamicList<point> recFPoints (recvFPoints);

            forAll(recFCells, rCell)
            {
                // get the cell label
                label recCell = recFCells[rCell];
                point recPoint = recFPoints[rCell];

                scalar uTau;
                if
				(
					uTauType_ == "interpPoint"
				 || uTauType_ == "effectiveDistance"
				)
                {
                    // Interpolate k
                    scalar kPoint = interpK->interpolate(recPoint, recCell);

                    // Compute friction velocity
                    uTau = Cmu25_*Foam::sqrt(kPoint);
                }
                else
                {
                    // Compute friction velocity
                    uTau = Cmu25_*Foam::sqrt(k[recCell]);
                }

                uTausToRetr[proci].append(uTau);
            }
        }
    }

    // Return
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendUTaus(proci, pBufsFCells);
            sendUTaus << uTausToRetr[proci];
        }
    }

    pBufsFCells.finishedSends();

    List<DynamicList<scalar>> uTausCmpl(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvUTaus(proci, pBufsFCells);
            DynamicList<scalar> recUTaus (recvUTaus);
            uTausCmpl[proci] = recUTaus;
        }
    }

    pBufsFCells.clear();

    // Complete uTau
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            forAll(uTausCmpl[proci], rCell)
            {
                scalar uTau = uTausCmpl[proci][rCell];
                label bLabel = bLabelsToRecv[proci][rCell];

                uTauAtIB_[Pstream::myProcNo()][bLabel] = uTau;
            }
        }
    }

    // Save
    saveUTau();
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::nutAtIB
(
    volScalarField& k,
    volScalarField& nu
)
{
    if (simulationType_ == "laminar")
    {
        // Loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            nutAtIB_[Pstream::myProcNo()][bCell] = 0.0;
        }
    }

    else if (nutWF_ == "nutkWallFunction")
    {
        // Loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Reset field
            nutAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // Get distance to the surface
            scalar yOrtho;
            if (useYEff_)
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            }

            // Get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // Compute yPlus
            scalar yPlus = uTau*yOrtho/nu[cellI];

            // Saves for later interpolation
            yPlusi_[cellI] = yPlus;

            // Compute the values at the surface
            if (yPlus > yPlusLam_)
            {
                nutAtIB_[Pstream::myProcNo()][bCell] =
					nu[cellI]*(yPlus*kappa_/Foam::log(E_*yPlus) - 1.0);
            }

            // Save
            nuti_[cellI] = nutAtIB_[Pstream::myProcNo()][bCell];
        }
    }
    else
    {
        FatalError
			<< nutWF_ << " condition for nut not implemented at the IB"
			<< exit(FatalError);
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::kAtIB
(
    List<scalar>& kIB,
    volScalarField& k,
    volScalarField& nu
)
{
    if (kWF_ == "kLowReWallFunction")
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // get distance to the surface
            scalar yOrtho;
            if (useYEff_)
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            }

            // Get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // Compute yPlus
            scalar yPlus = uTau*yOrtho/nu[cellI];

            // Saves for later interpolation
            yPlusi_[cellI] = yPlus;

            // Compute the values at the surface
            if (yPlus > yPlusLam_)
            {
                scalar Ck = -0.416;
                scalar Bk = 8.366;
                kIB[bCell] = (Ck/kappa_*Foam::log(yPlus) + Bk)*sqr(uTau);
            }
            else
            {
                scalar C = 11.0;
                scalar Cf = (1.0/sqr(yPlus + C) + 2.0*yPlus/pow3(C) - 1.0/sqr(C));
                kIB[bCell] = (2400.0/sqr(Ceps2_)*Cf)*sqr(uTau);
            }
        }

        // Ensure stability of computations
        forAll(kIB, bCell)
        {
            kIB[bCell] = max(kIB[bCell], small);
        }

        // Save
        forAll(kIB, bCell)
        {
            kAtIB_[Pstream::myProcNo()][bCell] = kIB[bCell];
        }
    }
    else
    {
        FatalError
			<< kWF_ << " condition for k not implemented at the IB"
			<< exit(FatalError);
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::omegaGAtIB
(
    List<scalar>& omegaIB,
    List<scalar>& GIB,
    volScalarField::Internal& G,
    const volVectorField& U,
    volScalarField& k,
    volScalarField& nu
)
{
    if (omegaWF_ == "omegaWallFunction")
    {
        // Blending switch
        bool blended = false; // should be an option

        // Load near wall dist
        //~ nearWallDist yWall(mesh_); // not used now

        // Get surface normal gradient
        List<DynamicList<vector>> snGradU;
        snGradU.setSize(Pstream::nProcs());
        snGradUAtIB(U, snGradU);

        // Loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Reset fields
            omegaIB[bCell] = 0.0;
            GIB[bCell] = 0.0;

            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // Get distance to the surface
            scalar yOrtho;
            if (useYEff_)
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            }

            // Compute magnitude of snGrad of U at the surface
            scalar magGradUWall = mag(snGradU[Pstream::myProcNo()][bCell]);

            // Get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // Compute local Reynolds number
            scalar Rey = yOrtho*uTau/nu[cellI];
            Rey /= Cmu25_;

            // Compute normalized variables
            const scalar yPlus = Cmu25_*Rey;
            const scalar uPlus = (1/kappa_)*Foam::log(E_*yPlus);

            // Saves for later interpolation
            yPlusi_[cellI] = yPlus;

            // Compute the values at the surface
            if (blended)
            {
                const scalar lamFrac = Foam::exp(-Rey/11);
                const scalar turbFrac = 1 - lamFrac;

                const scalar uStar = Foam::sqrt
                (
                    lamFrac*nu[cellI]*magGradUWall + turbFrac*sqr(uTau)
                );

                const scalar omegaVis = 6*nu[cellI]/(beta1_*Foam::sqr(yOrtho));
                const scalar omegaLog = uStar/(Cmu5_*kappa_*yOrtho);

                omegaIB[bCell] = lamFrac*omegaVis + turbFrac*omegaLog;
                GIB[bCell] =
				(
					lamFrac*G[cellI]
				  + turbFrac*sqr(uStar*magGradUWall*yOrtho/uPlus)
				   /(nu[cellI]*kappa_*yPlus)
				);
            }
            else
            {
                if (yPlus < yPlusLam_)
                {
                    omegaIB[bCell] = 6*nu[cellI]/(beta1_*Foam::sqr(yOrtho));
                    GIB[bCell] = G[cellI];
                }
                else
                {
                    const scalar uStar = uTau;

                    omegaIB[bCell] = uStar/(Cmu5_*kappa_*yOrtho);
                    GIB[bCell] =
					(
						sqr(uStar*magGradUWall*yOrtho/uPlus)
					   /(nu[cellI]*kappa_*yPlus)
					);
                }
            }
        }
    }
    else
    {
        FatalError
			<< omegaWF_
			<< " condition for omega and G not implemented at the IB"
			<< exit(FatalError);
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::epsilonGAtIB
(
    List<scalar>& epsilonIB,
    List<scalar>& GIB,
    volScalarField::Internal& G,
    const volVectorField& U,
    volScalarField& k,
    volScalarField& nu
)
{
    if (epsilonWF_ == "epsilonWallFunction")
    {
        // Load near wall dist
        //~ nearWallDist yWall(mesh_); // not used now

        // Get surface normal gradient
        List<DynamicList<vector>> snGradU;
        snGradU.setSize(Pstream::nProcs());
        snGradUAtIB(U, snGradU);

        // Loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // Reset fields
            epsilonIB[bCell] = 0.0;
            GIB[bCell] = 0.0;

            // Get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // Get distance to the surface
            scalar yOrtho;
            if (useYEff_)
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            }

            // Compute magnitude of snGrad of U at the surface
            scalar magGradUWall = mag(snGradU[Pstream::myProcNo()][bCell]);

            // Get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // Compute local Reynolds number
            scalar Rey = yOrtho*uTau/nu[cellI];
            Rey /= Cmu25_;

            // Compute normalized variables
            const scalar yPlus = Cmu25_*Rey;

            // Saves for later interpolation
            yPlusi_[cellI] = yPlus;

            if (yPlus > yPlusLam_)
            {
                epsilonIB[bCell] = pow3(uTau)/(kappa_*yOrtho);
                GIB[bCell] =
				(
					(nutAtIB_[Pstream::myProcNo()][bCell] + nu[cellI])
				   *magGradUWall*uTau/(kappa_*yOrtho)
				);
            }
            else
            {
                epsilonIB[bCell] = 2.0*k[cellI]*nu[cellI]/sqr(yOrtho);
                GIB[bCell] = G[cellI];
            }
        }
    }
    else
    {
        FatalError
			<< epsilonWF_
			<< " condition for epsilon and G not implemented at the IB"
			<< exit(FatalError);
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::calculateWallShearStress
(
    volVectorField& tauw,
    const volVectorField& U,
    const volScalarField& nu
)
{
    // Reset
    tauw *= 0.0;

    // Prepare list
    List<vector> tauwIB;
    tauwIB.setSize(boundaryCells_[Pstream::myProcNo()].size());

    // Prepare grad fields
    volTensorField gradU(fvc::grad(U));
    List<DynamicList<vector>> snGradU;
    snGradU.setSize(Pstream::nProcs());
    snGradUAtIB(U, snGradU);

    // Loop over boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // Get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // Calculate effective nu
        scalar nuEff = nu[cellI] + nutAtIB_[Pstream::myProcNo()][bCell];

        // Correct gradient
        vector normal = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;
        tensor correction =
		(
			normal
		   *(snGradU[Pstream::myProcNo()][bCell] - (normal & gradU[cellI]))
		);
        // Calculate dev tau
        symmTensor devTau = -nuEff*dev(twoSymm(gradU[cellI] + correction));

        // Calculate wall shear stress
        tauwIB[bCell] = -normal & devTau;

        // Save
        tauw[cellI] = tauwIB[bCell];
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::calculateForces
(
    volVectorField& fN,
    volVectorField& fT,
    volVectorField& tauw,
    const volScalarField& p,
    dictionary forceDict
)
{
    // Reset fields
    fN *= 0.0;
    fT *= 0.0;

    // Save for check
    volScalarField surfAdded
    (
        IOobject
        (
            "surfAdded",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );

    volVectorField surfPoints
    (
        IOobject
        (
            "surfPoints",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector::zero)
    );

    volVectorField tauws
    (
        IOobject
        (
            "tauws",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector::zero)
    );

    // Read dict
    scalar rhoInf = forceDict.lookupOrDefault<scalar>("rhoInf", 1000.0);
    scalar pRef = forceDict.lookupOrDefault<scalar>("pRef", 0.0);

    // Go through boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // Get cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

        // Save to added
		// Note (LK): type 0
        if (body_[outCellI] < 0.5 && body_[outCellI] >= thrSurf_)
        {
            surfAdded[outCellI] += 1.0;
        }
		// Note (LK): type 1
        else if (body_[outCellI] < thrSurf_ && body_[inCellI] < 1.0 - thrSurf_)
        {
            surfAdded[inCellI] += 1.0;
        }
		// Note (LK): type 2
        else
        {
            surfAdded[outCellI] += 1.0;
        }
    }

    // Go through boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // Get cell label
        label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
        label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

        // Decide where to put
        label whereI = outCellI;

        if (body_[outCellI] < thrSurf_ && body_[inCellI] < 1.0 - thrSurf_)
        {
            whereI = inCellI;
        }
        else
        {
            whereI = outCellI;
        }

        // Get and save surf point
        surfPoints[whereI] =
			boundaryCells_[Pstream::myProcNo()][bCell].sPoint_;

        // Add wall shear stress
        tauws[whereI] += tauw[outCellI]/surfAdded[whereI];

        // Get surface normal
        vector normal = -1*boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;

        // Get surface area
        scalar sA = boundaryCells_[Pstream::myProcNo()][bCell].sArea_;

        // calculate normal force
		// Note (LK): zero gradient considered
        fN[whereI] +=
			rhoInf*normal*sA*(p[outCellI] - pRef/rhoInf)/surfAdded[whereI];

        // calculate tangential force
		// Note (LK): minus in calculation of tauw
        fT[whereI] += -1*sA*rhoInf*tauw[outCellI]/surfAdded[whereI];
    }

    // Check if some surface cells were skipped
    forAll(surfaceCells_[Pstream::myProcNo()], sCell)
    {
        // Get cell label
        label cellI = surfaceCells_[Pstream::myProcNo()][sCell].sCell_;

        // Cell already added
        if (surfAdded[cellI] > 0.1)
        {
            continue;
        }

        // Prepare surf point and normal
        vector normal = -1*surfaceCells_[Pstream::myProcNo()][sCell].sNorm_;
        point surfPoint = surfaceCells_[Pstream::myProcNo()][sCell].sPoint_;

        // Get surface area
        scalar sA = surfaceCells_[Pstream::myProcNo()][sCell].sArea_;

        // Prepare total weight and value
        scalar totWeight = 0.0;
        vector totTauws = vector::zero;

        // Get values from neighbors
        forAll(mesh_.cells()[cellI], fI)
        {
            // Get face label
            label faceI = mesh_.cells()[cellI][fI];

            // Skip boundary faces
            if (!mesh_.isInternalFace(faceI))
            {
                continue;
            }

            // Get owner and neighbor
            label owner = mesh_.owner()[faceI];
            label neighbor = mesh_.neighbour()[faceI];

            // Get cell neighbor
            label nI(neighbor);
            if (neighbor == cellI)
            {
                nI = owner;
            }

            // Check if added from boundary cells
            if (surfAdded[nI] > 0.1)
            {
                point nSurfPoint = surfPoints[nI];
                scalar dist = mag(nSurfPoint - surfPoint);
                scalar weight = 1.0/dist;

                // Add
                totWeight += weight;
                totTauws += weight*tauws[nI];
            }
        }

        // Divide by total weight
        totTauws /= totWeight;

        // Calculate forces
        fN[cellI] += rhoInf*normal*sA*(p[cellI] - pRef/rhoInf);
        fT[cellI] += -1*sA*rhoInf*totTauws;

        // Check as added
        surfAdded[cellI] += 1.0;
    }

    // Note (LK): check what was added
    surfAdded.write();
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::calculateForceCoeffs
(
    scalar& Cl,
    scalar& Cd,
    volVectorField& fN,
    volVectorField& fT,
    dictionary forceDict
)
{
    // Read dict
    scalar rhoInf = forceDict.lookupOrDefault<scalar>("rhoInf", 1000.0);
    scalar magUInf = forceDict.lookupOrDefault<scalar>("magUInf", 1.0);
    scalar ARef = forceDict.lookupOrDefault<scalar>("ARef", 1.0);
    vector liftDir = forceDict.lookupOrDefault<vector>("liftDir", vector(0,1,0));
    vector dragDir = forceDict.lookupOrDefault<vector>("dragDir", vector(1,0,0));

    // Calculate dynamic pressure
    scalar pDyn = 0.5*rhoInf*magUInf*magUInf;

    // Calculate total force
    Field<vector> totForce = fN + fT;

    // Calculate coeffs fields
    Field<scalar> fieldCl((totForce & liftDir)/(ARef*pDyn));
    Field<scalar> fieldCd((totForce & dragDir)/(ARef*pDyn));

    // Calculate coefficients
    Cl = sum(fieldCl);
    Cd = sum(fieldCd);
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::saveUTau()
{
    // Reset saved data
    uTaui_ *= 0.0;

    // Loop over boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // Get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // Save
        uTaui_[cellI] = uTauAtIB_[Pstream::myProcNo()][bCell];
    }
}

// ------------------------------------------------------------------------- //

void ibDirichletBCs::snGradUAtIB
(
    const volVectorField& U,
    List<DynamicList<vector>>& snGradU
)
{
    // Note (LK): possibility to add new sn grad schemes here

    // Loop over boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // Get the cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // Get distance to the surface
        scalar yOrtho;
        if (useYEff_)
        {
            yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
        }
        else
        {
            yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
        }

        // Calculate surface normal gradient
		// Note (LK): not moving solid considered, should be changed
        vector snGrad = (vector::zero - U[cellI])/yOrtho;

        // Assign
        snGradU[Pstream::myProcNo()].append(snGrad);
    }
}


// ************************************************************************* //
