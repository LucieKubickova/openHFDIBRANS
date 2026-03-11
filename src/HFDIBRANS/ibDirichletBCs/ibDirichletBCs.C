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

#include "ibDirichletBCs.H"

using namespace Foam;

//---------------------------------------------------------------------------//
ibDirichletBCs::ibDirichletBCs
(
    const fvMesh& mesh,
    const volScalarField& body,
    List<DynamicList<boundaryCell>>& boundaryCells,
    labelField& isBoundaryCell
)
:
mesh_(mesh),
body_(body),
boundaryCells_(boundaryCells),
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
transportProperties_
(
    IOobject
    (
        "transportProperties",
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
tauwi_
(
    IOobject
    (
        "tauwi",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedVector("zero", dimless, vector::zero)
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
Prt_(0.85),
Cmu_(0.09),
Ceps2_(1.9),
beta1_(0.075)
{
    // initiate lists
    nutAtIB_.setSize(Pstream::nProcs());
    kAtIB_.setSize(Pstream::nProcs());
    uTauAtIB_.setSize(Pstream::nProcs());

    // read turbulence properties
    turbulenceProperties_.lookup("simulationType") >> simulationType_;
    
    // read HFDIBDEMDict
    thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));
    useYEff_ = HFDIBDEMDict_.lookupOrDefault<bool>("useEffectiveDist", true);
    uTauType_ = HFDIBDEMDict_.lookupOrDefault<word>("uTauType", "freeStreamCell");
    uTauCoeff_ = HFDIBDEMDict_.lookupOrDefault<scalar>("uTauCoeff", 1.0);

    // read boundary condition for velocity
    HFDIBBCsDict_ = HFDIBDEMDict_.subDict("wallFunctions");
    UBC_ = HFDIBBCsDict_.lookupOrDefault<word>("U", "noSlip");

    // read simulation type
    if (simulationType_ != "laminar")
    {
        nutWF_ = HFDIBBCsDict_.lookupOrDefault<word>("nut", "undefined");
        alphatWF_ = HFDIBBCsDict_.lookupOrDefault<word>("alphat", "undefined");
        kWF_ = HFDIBBCsDict_.lookupOrDefault<word>("k", "undefined");
        omegaWF_ = HFDIBBCsDict_.lookupOrDefault<word>("omega", "undefined");
        epsilonWF_ = HFDIBBCsDict_.lookupOrDefault<word>("epsilon", "undefined");
    }

    // compute turbulence parameters
    Cmu75_ = Foam::pow(Cmu_, 0.75);
    Cmu25_ = pow025(Cmu_);
    Cmu5_ = Foam::sqrt(Cmu_);
    calcYPlusLam();
}

//---------------------------------------------------------------------------//
ibDirichletBCs::~ibDirichletBCs()
{
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::calcYPlusLam
(
)
{
    yPlusLam_ = 11.0;

    for (int i=0; i<10; i++)
    {
        yPlusLam_ = Foam::log(max(E_*yPlusLam_, 1))/kappa_;
    }
}

//---------------------------------------------------------------------------//
scalar ibDirichletBCs::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11;
    scalar tolerance = 0.01;
    label maxIters = 10;

    for (int iter = 0; iter < maxIters; ++iter)
    {
        const scalar f = ypt - (Foam::log(E_*ypt)/kappa_ + P)/Prat;
        const scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
        const scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::setSizeToLists
(
)
{
    if (simulationType_ != "laminar")
    {
        // set size
        nutAtIB_[Pstream::myProcNo()].setSize(boundaryCells_[Pstream::myProcNo()].size());
        kAtIB_[Pstream::myProcNo()].setSize(boundaryCells_[Pstream::myProcNo()].size());
        uTauAtIB_[Pstream::myProcNo()].setSize(boundaryCells_[Pstream::myProcNo()].size());

        // reset
        nutAtIB_[Pstream::myProcNo()] = 0.0;
        kAtIB_[Pstream::myProcNo()] = 0.0;
        uTauAtIB_[Pstream::myProcNo()] = 0.0;
    }
}

//---------------------------------------------------------------------------//
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
            // assign zero
            UIB[uCell] = ibZero(UIB[uCell]);
        }
    }

    else if (UBC_ == "partialSlip")
    {
        // read partial clip coefficient
        scalar alpha = HFDIBBCsDict_.lookupOrDefault<scalar>("UCoeff", 0.0);

        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // reset field
            UIB[bCell] = vector::zero;

            // get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // get surface normal
            vector normal = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;

            // calculate value at boundary
            UIB[bCell] += (1 - alpha)*transform(I - sqr(normal), U[cellI]);
        }
    }

    else
    {
        FatalError << UBC_ << " condition for U in " << simulationType_ << " not implemented at the IB" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::TAtIB
(
    List<scalar>& TIB,
    scalar TIn
)
{
    forAll(TIB, bCell)
    {
        // assign value
        TIB[bCell] = TIn;
    }
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::updateUTauAtIB
(
    volScalarField& k
)
{
    // prepare interpolation scheme
    dictionary HFDIBInnerSchemes = fvSchemes_.subDict("HFDIBSchemes").subDict("innerSchemes");
    autoPtr<interpolation<scalar>> interpK = interpolation<scalar>::New(HFDIBInnerSchemes, k);

    // prepare sync
    List<DynamicList<label>> fCellsToSync(Pstream::nProcs());
    List<DynamicList<point>> fPointsToSync(Pstream::nProcs());
    List<DynamicList<label>> bLabelsToRecv(Pstream::nProcs());

    // if uTau from boundary cell
    if (uTauType_ == "boundaryCell")
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // reset field
            uTauAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // compute uTau
            uTauAtIB_[Pstream::myProcNo()][bCell] = Cmu25_*Foam::sqrt(k[cellI]);
        }
    }

    else if (uTauType_ == "effectiveDistance" or uTauType_ == "cellDistance" or uTauType_ == "coeffDistance")
    {
        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // reset field
            uTauAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
            label fCell1 = boundaryCells_[Pstream::myProcNo()][bCell].fCell1_;
            label fCell2 = boundaryCells_[Pstream::myProcNo()][bCell].fCell2_;
            label fProc1 = boundaryCells_[Pstream::myProcNo()][bCell].fProc1_;
            label fProc2 = boundaryCells_[Pstream::myProcNo()][bCell].fProc2_;

            // get surface point and normal
            point sPoint = boundaryCells_[Pstream::myProcNo()][bCell].sPoint_;
            vector sNorm = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;

            // get distance
            scalar dist(0.0);
            if (uTauType_ == "effectiveDistance")
            {
                dist = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else if (uTauType_ == "cellDistance")
            {
                dist = Foam::pow(mesh_.V()[cellI],0.333);
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
            if (pointInCell(distPoint, cellI))
            {
                kCell = cellI;
                kProc = Pstream::myProcNo();
            }
            else //~ if (pointInCell(yEffPoint, fCell1)) // Note (LK): the distPoint should not be farther than one cell away
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
                // interpolate k
                scalar kPoint = interpK->interpolate(distPoint, kCell);

                // compute friction velocity
                uTauAtIB_[Pstream::myProcNo()][bCell] = Cmu25_*Foam::sqrt(kPoint);
            }

            else
            {
                fCellsToSync[kProc].append(kCell);
                fPointsToSync[kProc].append(distPoint);
                bLabelsToRecv[kProc].append(bCell);
            }
        }
    }

    else if (uTauType_ == "freeStreamCell" or uTauType_ == "interpPoint")
    {
        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // reset field
            uTauAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // prepare
            label fCell(0);
            label fProc(-1);
            point fPoint(vector::zero);

            // use the first interpolation point
            fCell = boundaryCells_[Pstream::myProcNo()][bCell].fCell1_;
            fProc = boundaryCells_[Pstream::myProcNo()][bCell].fProc1_;
            fPoint = boundaryCells_[Pstream::myProcNo()][bCell].fPoint1_;

            // check body -- Note (LK): not really working for pipe or anything, but idea to smooth out
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

            // compute uTau based on values from the free stream
            if (Pstream::myProcNo() == fProc)
            {
                if (uTauType_ == "interpPoint")
                {
                    // interpolate k
                    scalar kPoint = interpK->interpolate(fPoint, fCell);

                    // compute friction velocity
                    uTauAtIB_[Pstream::myProcNo()][bCell] = Cmu25_*Foam::sqrt(kPoint);
                }

                else
                {
                    // compute friction velocity
                    uTauAtIB_[Pstream::myProcNo()][bCell] = Cmu25_*Foam::sqrt(k[fCell]);
                }
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
        FatalError << "uTau calculation type " << uTauType_ << " not implemented" << exit(FatalError);
    }

    // sync with other processors
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

    // recieve
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
                if (uTauType_ == "firstInterpPoint" or uTauType_ == "effectiveDistance")
                {
                    // interpolate k
                    scalar kPoint = interpK->interpolate(recPoint, recCell);

                    // compute friction velocity
                    uTau = Cmu25_*Foam::sqrt(kPoint);
                }

                else
                {
                    // compute friction velocity
                    uTau = Cmu25_*Foam::sqrt(k[recCell]);
                }

                uTausToRetr[proci].append(uTau);
            }
        }
    }

    // return
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

    // complete uTau
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

    // post processing
    postProcessUTau();

    // save
    saveUTau();
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::nutAtIB
(
    volScalarField& k,
    volScalarField& nu
)
{
    if (nutWF_ == "nutkWallFunction")
    {
        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // reset field
            nutAtIB_[Pstream::myProcNo()][bCell] = 0.0;

            // get cell label
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

            // get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // compute yPlus
            scalar yPlus = uTau*yOrtho/nu[cellI];

            // saves for later interpolation
            yPlusi_[cellI] = yPlus;

            // compute the values at the surface
            if (yPlus > yPlusLam_)
            {
                nutAtIB_[Pstream::myProcNo()][bCell] = nu[cellI]*(yPlus*kappa_/Foam::log(E_*yPlus) - 1.0);
            }

            // save
            nuti_[cellI] = nutAtIB_[Pstream::myProcNo()][bCell];
        }
    }

    else
    {
        FatalError << nutWF_ << " condition for nut not implemented at the IB" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::correctAlphatAtIB
(
    List<scalar>& alphatIB,
    const volScalarField& nu
)
{
    if (alphatWF_ == "alphatJayatillekeWallFunction")
    {
        // Molecular Prandtl number
        const scalar Pr
        (
            dimensionedScalar("Pr", dimless, transportProperties_).value()
        );

        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get cell label
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

            // get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // compute yPlus
            scalar yPlus = uTau*yOrtho/nu[cellI];

            // saves for later interpolation
            yPlusi_[cellI] = yPlus;

            // Molecular-to-turbulent Prandtl number ratio
            const scalar Prat = Pr/Prt_;

            // Thermal sublayer thickness
            const scalar P = 9.24*(Foam::pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*Foam::exp(-0.007*Prat)); // Psmooth function in source
            const scalar yPlusTherm = this->yPlusTherm(P, Prat);

            // Update turbulent thermal conductivity
            if (yPlus > yPlusTherm)
            {
                const scalar kt =
                    nu[cellI]*(yPlus/(Prt_*(Foam::log(E_*yPlus)/kappa_ + P)) - 1.0/Pr);

                alphatIB[bCell] = max(scalar(0), kt);
            }
            else
            {
                alphatIB[bCell] = 0.0;
            }
        }
    }

    else if (simulationType_ == "laminar")
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            alphatIB[bCell] = 0.0;
        }
    }
}

//---------------------------------------------------------------------------//
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
            // get cell label
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

            // get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // compute yPlus
            scalar yPlus = uTau*yOrtho/nu[cellI];

            // saves for later interpolation
            yPlusi_[cellI] = yPlus;

            // compute the values at the surface
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

        // ensure stability of computations
        forAll(kIB, bCell)
        {
            kIB[bCell] = max(kIB[bCell], SMALL);
        }

        // save
        forAll(kIB, bCell)
        {
            kAtIB_[Pstream::myProcNo()][bCell] = kIB[bCell];
        }

        // save kIB
        //~ word fileName = "k.dat";
        //~ word outDir = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/ZZ_python";

        //~ // prepare outFile
        //~ autoPtr<OFstream> outFilePtr;
        //~ outFilePtr.reset(new OFstream(outDir/fileName));
        //~ outFilePtr() << "cellI,x,y,z,V,k" << endl;

        //~ // loop over cells
        //~ forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        //~ {
            //~ // get cell label
            //~ label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            //~ // get distance to the surface
            //~ scalar yOrtho;
            //~ if (useYEff_)
            //~ {
            //~     yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            //~ }
            //~ else
            //~ {
            //~     yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            //~ }

            //~ // get coordinates and volume
            //~ scalar x = mesh_.C()[cellI].x();
            //~ scalar y = mesh_.C()[cellI].y();
            //~ if (y < 0.0) // UGLYYYYYYYYYYYY
            //~ {
                //~ y -= yOrtho;
            //~ }
            //~ else
            //~ {
                //~ y += yOrtho;
            //~ }
            //~ scalar z = mesh_.C()[cellI].z();
            //~ scalar V = mesh_.V()[cellI];

            //~ // get the fields value
            //~ scalar kk = kIB[bCell];

            //~ // write
            //~ outFilePtr() << cellI
                //~ << "," << x
                //~ << "," << y
                //~ << "," << z
                //~ << "," << V
                //~ << "," << kk
                //~ << endl;
        //~ }
    }

    else
    {
        FatalError << kWF_ << " condition for k not implemented at the IB" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
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
        // blending switch
        bool blended(false); // should be an option

        // load near wall dist
        //~ nearWallDist yWall(mesh_); // not used now

        // get surface normal gradient
        List<DynamicList<vector>> snGradU;
        snGradU.setSize(Pstream::nProcs());
        snGradUAtIB(U, snGradU);

        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // reset fields
            omegaIB[bCell] = 0.0;
            GIB[bCell] = 0.0;

            // get cell label
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

            // compute magnitude of snGrad of U at the surface
            scalar magGradUWall = mag(snGradU[Pstream::myProcNo()][bCell]);

            // get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // compute local Reynolds number
            scalar Rey = yOrtho*uTau/nu[cellI];
            Rey /= Cmu25_;
            
            // compute normalized variables
            const scalar yPlus = Cmu25_*Rey;
            const scalar uPlus = (1/kappa_)*Foam::log(E_*yPlus);

            // saves for later interpolation
            yPlusi_[cellI] = yPlus;

            // compute the values at the surface
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
                GIB[bCell] = lamFrac*G[cellI] + turbFrac*sqr(uStar*magGradUWall*yOrtho/uPlus)/(nu[cellI]*kappa_*yPlus);
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
                    GIB[bCell] = sqr(uStar*magGradUWall*yOrtho/uPlus)/(nu[cellI]*kappa_*yPlus);
                }
            }
        }
    }

    else
    {
        FatalError << omegaWF_ << " condition for omega and G not implemented at the IB" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
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
        // load near wall dist
        //~ nearWallDist yWall(mesh_); // not used now

        // get surface normal gradient
        List<DynamicList<vector>> snGradU;
        snGradU.setSize(Pstream::nProcs());
        snGradUAtIB(U, snGradU);

        // loop over boundary cells
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // reset fields
            epsilonIB[bCell] = 0.0;
            GIB[bCell] = 0.0;

            // get cell label
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

            // compute magnitude of snGrad of U at the surface
            scalar magGradUWall = mag(snGradU[Pstream::myProcNo()][bCell]);

            // get the friction velocity
            scalar uTau = uTauAtIB_[Pstream::myProcNo()][bCell];

            // compute local Reynolds number
            scalar Rey = yOrtho*uTau/nu[cellI];
            Rey /= Cmu25_;

            // compute normalized variables
            const scalar yPlus = Cmu25_*Rey;

            // saves for later interpolation
            yPlusi_[cellI] = yPlus;

            if (yPlus > yPlusLam_)
            {
                epsilonIB[bCell] = pow3(uTau)/(kappa_*yOrtho);
                GIB[bCell] = (nutAtIB_[Pstream::myProcNo()][bCell] + nu[cellI])*magGradUWall*uTau/(kappa_*yOrtho);
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
        FatalError << epsilonWF_ << " condition for epsilon and G not implemented at the IB" << exit(FatalError);
    }
}


//---------------------------------------------------------------------------//
// Note (LK): not parallelized, but not really used either
void ibDirichletBCs::postProcessUTau
(
)
{
    //~ // initialize
    //~ scalar totA(0.0);

    //~ // loop over faces
    //~ forAll(mesh_.cells()[cellI], fI)
    //~ {
        //~ // get face label
        //~ label faceI = mesh_.cells()[cellI][fI];

        //~ // skip boundary faces
        //~ if (faceI >= mesh_.faceNeighbour().size())
        //~ {
            //~ continue;
        //~ }

        //~ // get cell labels
        //~ label owner = mesh_.faceOwner()[faceI];
        //~ label neighbor = mesh_.faceNeighbour()[faceI];

        //~ // skip in-solid cells
        //~ if (body_[owner] >= 0.5 or body_[neighbor] >= 0.5)
        //~ {
            //~ continue;
        //~ }

        //~ // get uTau values
        //~ scalar uTauO = Cmu25_*Foam::sqrt(k[owner]);
        //~ scalar uTauN = Cmu25_*Foam::sqrt(k[neighbor]);

        //~ // calculate the average value
        //~ uTauAtIB_[Pstream::myProcNo()][bCell] += mag(mesh_.Sf()[faceI])*(uTauO*mag(mesh_.Cf()[faceI]-mesh_.C()[neighbor]) + uTauN*mag(mesh_.Cf()[faceI]-mesh_.C()[owner]))/(mag(mesh_.Cf()[faceI] - mesh_.C()[neighbor]) + mag(mesh_.Cf()[faceI] - mesh_.C()[owner]));

        //~ // add face area to total
        //~ totA += mag(mesh_.Sf()[faceI]);
    //~ }

    //~ // loop over boundary cells -- does not make a difference
    //~ forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    //~ {
        //~ // get the cell label
        //~ label cellN = boundaryCells_[Pstream::myProcNo()][nCell].bCell_;

        //~ // flag
        //~ bool isNeigh = false;

        //~ // loop over neighbor cells
        //~ forAll(mesh_.cellCells()[cellN], cN)
        //~ {
            //~ if (cellI == mesh_.cellCells()[cellN][cN])
            //~ {
                //~ isNeigh = true;
                //~ break;
            //~ }
        //~ }

        //~ // for neighbors add to uTau
        //~ if (isNeigh)
        //~ {
            //~ // find common face
            //~ label faceI;
            //~ forAll(mesh_.cells()[cellI], fI)
            //~ {
                //~ // flag
                //~ bool found = false;

                //~ // get the face label
                //~ faceI = mesh_.cells()[cellI][fI];

                //~ forAll(mesh_.cells()[cellN], fN)
                //~ {
                    //~ // get the face label
                    //~ label faceN = mesh_.cells()[cellN][fN];

                    //~ if (faceI == faceN)
                    //~ {
                        //~ found = true;
                        //~ break;
                    //~ }
                //~ }

                //~ if (found)
                //~ {
                    //~ break;
                //~ }
            //~ }

            //~ // compute uTau for each cell
            //~ scalar uTauO = Cmu25_*Foam::sqrt(k[cellI]);
            //~ scalar uTauN = Cmu25_*Foam::sqrt(k[cellN]);

            //~ // add to total uTau
            //~ uTauAtIB_[Pstream::myProcNo()][bCell] += mag(mesh_.Sf()[faceI])*(uTauO*mag(mesh_.Cf()[faceI]-mesh_.C()[cellN]) + uTauN*mag(mesh_.Cf()[faceI]-mesh_.C()[cellI]))/(mag(mesh_.Cf()[faceI] - mesh_.C()[cellN]) + mag(mesh_.Cf()[faceI] - mesh_.C()[cellI]));

            //~ // add face area to total
            //~ totA += mag(mesh_.Sf()[faceI]);
        //~ }
    //~ }

    //~ // divide by total area
    //~ uTauAtIB_[Pstream::myProcNo()][bCell] /= totA;
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::calculateWallShearStress
(
    const volVectorField& U,
    const volScalarField& nu
)
{
    // reset
    tauwi_ *= 0.0;

    // prepare list
    List<vector> tauwIB;
    tauwIB.setSize(boundaryCells_[Pstream::myProcNo()].size());

    // prepare grad fields
    volTensorField gradU = fvc::grad(U);
    List<DynamicList<vector>> snGradU;
    snGradU.setSize(Pstream::nProcs());
    snGradUAtIB(U, snGradU);

    // loop over boundar cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // calculate effective nu
        scalar nuEff = nu[cellI] + nutAtIB_[Pstream::myProcNo()][bCell];

        // correct gradient
        vector normal = boundaryCells_[Pstream::myProcNo()][bCell].sNorm_;
        tensor correction = normal * (snGradU[Pstream::myProcNo()][bCell] - (normal & gradU[cellI]));

        // calculate dev tau
        symmTensor devTau = -nuEff * dev(twoSymm(gradU[cellI] + correction));

        // calculate wall shear stress
        tauwIB[bCell] = -normal & devTau;

        // save
        tauwi_[cellI] = tauwIB[bCell];
    }
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::saveUTau
(
)
{
    // reset saved data
    uTaui_ *= 0.0;

    // loop over boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // save
        uTaui_[cellI] = uTauAtIB_[Pstream::myProcNo()][bCell];
    }
}

//---------------------------------------------------------------------------//
bool ibDirichletBCs::pointInCell // copy from lineIntInfo
(
    point pToCheck,
    label cToCheck
)
{
    const labelList& cellFaces(mesh_.cells()[cToCheck]);
    forAll(cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        vector outNorm = mesh_.Sf()[fI];
        outNorm = (mesh_.faceOwner()[fI] == cToCheck) ? outNorm : (-1*outNorm);

        if (((pToCheck - mesh_.Cf()[fI]) & outNorm) > 0)
        {
            return false;
        }
    }
    return true;
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::snGradUAtIB
(
    const volVectorField& U,
    List<DynamicList<vector>>& snGradU
)
{
    // Note (LK): possibility to add new sn grad schemes here

    // loop over boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get the cell label
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

        // calculate surface normal gradient
        vector snGrad = (vector::zero - U[cellI])/yOrtho; // Note (LK): not moving solid considered, should be changed

        // assign
        snGradU[Pstream::myProcNo()].append(snGrad);
    }
}

// ************************************************************************* //
