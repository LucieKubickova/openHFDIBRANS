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
#include "openHFDIBRANS.H"

using namespace Foam;

//---------------------------------------------------------------------------//
openHFDIBRANS::openHFDIBRANS
(
    const fvMesh& mesh,
    const volScalarField& body
) :
mesh_(mesh),
body_(body),
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
)
{
    // initiate lists
    boundaryCells_.setSize(Pstream::nProcs());
    surfaceCells_.setSize(Pstream::nProcs());
    internalCells_.setSize(Pstream::nProcs());

    // initialize classes
    ibInterpolation_.set(new ibInterpolation(mesh_, body_, boundaryCells_, surfaceCells_, internalCells_, isBoundaryCell_));
    ibDirichletBCs_.set(new ibDirichletBCs(mesh_, body_, boundaryCells_, isBoundaryCell_));

    // read HFDIBDEM dictionary
    save_ = HFDIBDEMDict_.lookupOrDefault<bool>("saveIntInfo", false);
    cpDisToInner_ = HFDIBDEMDict_.lookupOrDefault<bool>("copyDisToInner", false);
    scaleDisG_ = HFDIBDEMDict_.lookupOrDefault<bool>("scaleDisG", false);
    scaleCoeff_ = HFDIBDEMDict_.lookupOrDefault<scalar>("scaleCoeff", 1.0);
    useYEff_ = HFDIBDEMDict_.lookupOrDefault<bool>("useEffectiveDist", false);
    thrSurf_ = readScalar(HFDIBDEMDict_.lookup("surfaceThreshold"));

    // read fvSchemes
    HFDIBOuterSchemes_ = fvSchemes_.subDict("HFDIBSchemes").subDict("outerSchemes");

    // identify boundary cells
    ibInterpolation_->findBoundaryCells();
    ibInterpolation_->findSurfaceCells();

    // set size to lists
    ibDirichletBCs_->setSizeToLists();

    // compute distance to the immersed boundary
    ibInterpolation_->calculateBoundaryDist();
    ibInterpolation_->calculateSurfaceDist();

    // calculate interpolation points
    ibInterpolation_->calculateInterpolationPoints();

    // save data
    if (save_)
    {
        // save boundary cells as cell sets
        ibInterpolation_->saveBoundaryCells();
        ibInterpolation_->saveSurfaceCells();

        // create output directory to save data for python
        word outDir = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/ZZ_python";
        if (!isDir(outDir))
        {
            mkDir(outDir);
        }

        // save interpolation data
        ibInterpolation_->saveInterpolationInfo(outDir, "interpolationInfo");
    }
}

//---------------------------------------------------------------------------//
openHFDIBRANS::~openHFDIBRANS()
{
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::computeUi
(
    volVectorField& U,
    volVectorField& Ui,
    word surfType // NOTE: add UIn
)
{
    // reset imposed field
    Ui *= 0.0;

    // calculate values at the immersed boundary
    List<vector> UIB;
    if (surfType == "lambdaBased")
    {
        UIB.setSize(surfaceCells_[Pstream::myProcNo()].size());
    }
    else
    {
        UIB.setSize(boundaryCells_[Pstream::myProcNo()].size());
    }

    ibDirichletBCs_->UAtIB(UIB, "noSlip");

    // get references
    volScalarField& yPlusi = ibDirichletBCs_->getYPlusi();
    scalar yPlusLam = ibDirichletBCs_->getYPlusLam();
        
    // calculate log scales for interpolation
    List<scalar> logScales;
    logScales.setSize(UIB.size());

    if (surfType != "lambdaBased")
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // get distance to surface
            scalar yOrtho;
            if (useYEff_)
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            }

            // calculate the local log scale
            logScales[bCell] = yPlusi[cellI]/yOrtho*ibDirichletBCs_->getE();
        }
    }

    // read interpolation schemes from fvSchemes
    ITstream UIBScheme = HFDIBOuterSchemes_.lookup("U");
    word interpType = UIBScheme[0].wordToken();

    // interpolation
    if (interpType == "unifunctional")
    {
        ibInterpolation_->unifunctionalInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales);
    }

    else if (interpType == "lambdaBased")
    {
        ibInterpolation_->lambdaBasedInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales);
    }

    else if (interpType == "switched")
    {
        ibInterpolation_->switchedInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "outerInner")
    {
        ibInterpolation_->outerInnerInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "inner")
    {
        ibInterpolation_->innerInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales, yPlusi, yPlusLam);
    }

    else
    {
        FatalError << "Interpolation type " << UIBScheme << " for field U not implemented" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::computeKi
(
    volScalarField& k,
    volScalarField& ki,
    volScalarField& nu
)
{
    // prepare lists
    List<scalar> kIB;
    kIB.setSize(boundaryCells_[Pstream::myProcNo()].size());

    // compute values at the immersed boundary
    ibDirichletBCs_->kAtIB(kIB, k, nu);

    // get references
    volScalarField& yPlusi = ibDirichletBCs_->getYPlusi();
    scalar yPlusLam = ibDirichletBCs_->getYPlusLam();

    // calculate log scales for interpolation
    List<scalar> logScales;
    logScales.setSize(boundaryCells_[Pstream::myProcNo()].size());

    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // get distance to surface
        scalar yOrtho;
        if (useYEff_)
        {
            yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
        }
        else
        {
            yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
        }

        // calculate the local log scale
        logScales[bCell] = yPlusi[cellI]/yOrtho;
    }

    // read interpolation schemes from fvSchemes
    ITstream kIBScheme = HFDIBOuterSchemes_.lookup("k");
    word interpType = kIBScheme[0].wordToken();

    // interpolation
    if (interpType == "unifunctional")
    {
        ibInterpolation_->unifunctionalInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales);
    }

    else if (interpType == "switched")
    {
        ibInterpolation_->switchedInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "outerInner")
    {
        ibInterpolation_->outerInnerInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "inner")
    {
        ibInterpolation_->innerInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales, yPlusi, yPlusLam);
    }

    else
    {
        FatalError << "Interpolation type " << kIBScheme << " for field k not implemented" << exit(FatalError);
    }

    // bound ki
    forAll(ki, cellI)
    {
        ki[cellI] = max(ki[cellI], small);
    }

    // TODO: blended interpolation
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::computeTi
(
    volScalarField& T,
    volScalarField& Ti,
    word surfType,
    scalar TIn
)
{
    // reset imposed field
    Ti *= 0.0;

    // assign the values in in-solid cells
    forAll(body_, cellI)
    {
        if (body_[cellI] >= 0.5)
        {
            Ti[cellI] = TIn;
        }
    }

    // calculate values at the immersed boundary
    List<scalar> TIB;
    if (surfType == "lambdaBased")
    {
        TIB.setSize(surfaceCells_[Pstream::myProcNo()].size());
    }
    else
    {
        TIB.setSize(boundaryCells_[Pstream::myProcNo()].size());
    }

    // compute values at the immersed boundary
    ibDirichletBCs_->TAtIB(TIB, TIn);

    // get references
    volScalarField& yPlusi = ibDirichletBCs_->getYPlusi();
    scalar yPlusLam = ibDirichletBCs_->getYPlusLam();
        
    // calculate log scales for interpolation
    List<scalar> logScales;
    logScales.setSize(TIB.size());

    if (surfType != "lambdaBased")
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get cell label
            label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // get distance to surface
            scalar yOrtho;
            if (useYEff_)
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            }
            else
            {
                yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            }

            // calculate the local log scale
            logScales[bCell] = yPlusi[cellI]/yOrtho*ibDirichletBCs_->getE();
        }
    }

    // read interpolation schemes from fvSchemes
    ITstream TIBScheme = HFDIBOuterSchemes_.lookup("T");
    word interpType = TIBScheme[0].wordToken();

    // use boundary condition
    if (interpType == "unifunctional")
    {
        ibInterpolation_->unifunctionalInterp<scalar, volScalarField>(TIBScheme, T, Ti, TIB, logScales);
    }

    else if (interpType == "lambdaBased")
    {
        ibInterpolation_->lambdaBasedInterp<scalar, volScalarField>(TIBScheme, T, Ti, TIB, logScales);
    }

    else if (interpType == "switched")
    {
        ibInterpolation_->switchedInterp<scalar, volScalarField>(TIBScheme, T, Ti, TIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "outerInner")
    {
        ibInterpolation_->outerInnerInterp<scalar, volScalarField>(TIBScheme, T, Ti, TIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "inner")
    {
        ibInterpolation_->innerInterp<scalar, volScalarField>(TIBScheme, T, Ti, TIB, logScales, yPlusi, yPlusLam);
    }

    else
    {
        FatalError << "Interpolation type " << TIBScheme << " for field T not implemented" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::updateUTau
(
    volScalarField& k
)
{
    ibDirichletBCs_->updateUTauAtIB(k);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctNut
(
    volScalarField& k,
    volScalarField& nu
)
{
    ibDirichletBCs_->correctNutAtIB(k, nu);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctOmegaG
(
    volScalarField& omega,
    volScalarField::Internal& G,
    const volVectorField& U,
    volScalarField& k,
    volScalarField& nu,
    volScalarField& surface
)
{
    // prepare lists
    List<scalar> omegaIB;
    List<scalar> GIB;

    omegaIB.setSize(boundaryCells_[Pstream::myProcNo()].size());
    GIB.setSize(boundaryCells_[Pstream::myProcNo()].size());

    // calculate values at the immersed boundary
    ibDirichletBCs_->omegaGAtIB(omegaIB, GIB, G, U, k, nu);

    // omega scaling
    if (scaleDisG_)
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get cell label
            //~ label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // get cell scales
            //~ scalar yOrtho;
            //~ if (useYEff_)
            //~ {
                //~ yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            //~ }
            //~ else
            //~ {
                //~ yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            //~ }
            //~ scalar V = mesh_.V()[cellI];
            //~ scalar l = Foam::pow(V, 0.333);

            // user-defined scaling
            omegaIB[bCell] = omegaIB[bCell]*scaleCoeff_;
            GIB[bCell] = GIB[bCell]*scaleCoeff_;

            // old versions
            //~ omegaIB[bCell] = omegaIB[bCell]*(2*yOrtho)/l;
            //~ GIB[bCell] = GIB[bCell]*(2*yOrtho)/l;
            //~ omegaIB[bCell] = omegaIB[bCell]*l/(2*yOrtho); // scalingRight
            //~ GIB[bCell] = GIB[bCell]*l/(2*yOrtho);

            // assign
            //~ if (body_[cellI] > thrSurf_) // yOrtho < l
            //~ {
                //~ omegaIB[bCell] = omegaIB[bCell]*l/(2*yOrtho); // scalingRight
                //~ GIB[bCell] = GIB[bCell]*l/(2*yOrtho);
            //~ }

            //~ else // yOrtho > l
            //~ {
                //~ omegaIB[bCell] = omegaIB[bCell]/l*(2*yOrtho); // scalingRightInverse
                //~ GIB[bCell] = GIB[bCell]/l*(2*yOrtho);
                //~ //~ omegaIB[bCell] = omegaIB[bCell]/l*(yOrtho); // NOTE: Martins function
                //~ //~ GIB[bCell] = GIB[bCell]/l*(yOrtho);
            //~ }
        }
    }

    // assign the values for boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // assign
        omega[cellI] = omegaIB[bCell];
        G[cellI] = GIB[bCell];
    }

    // sync boundary conditions
    omega.correctBoundaryConditions();
    //~ G.correctBoundaryConditions();

    // calculate maximum omega
    scalar inOmega = 0.0;
    inOmega = max(omegaIB); // internal patch fields for walls should be included as well
    reduce(inOmega, maxOp<scalar>());

    // assign the values in in-solid cells
    forAll(surface, cellI)
    {
        if (surface[cellI] == 1.0)
        {
            if (body_[cellI] >= 0.5)
            {
                omega[cellI] = inOmega;
                G[cellI] = 0.0;
            }
        }
    }

    // sync boundary conditions
    omega.correctBoundaryConditions();

    // correct the inner boundary cells
    if (cpDisToInner_)
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get cell labels
            label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
            label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

            // assign
            omega[inCellI] = omega[outCellI];
            G[inCellI] = G[outCellI];
        }
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctEpsilonG
(
    volScalarField& epsilon,
    volScalarField::Internal& G,
    const volVectorField& U,
    volScalarField& k,
    volScalarField& nu,
    volScalarField& surface
)
{
    // prepare lists
    List<scalar> epsilonIB;
    List<scalar> GIB;

    epsilonIB.setSize(boundaryCells_[Pstream::myProcNo()].size());
    GIB.setSize(boundaryCells_[Pstream::myProcNo()].size());

    // calculate values at the immersed boundary
    ibDirichletBCs_->epsilonGAtIB(epsilonIB, GIB, G, U, k, nu);

    // epsilon scaling
    if (scaleDisG_)
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get cell label
            //~ label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

            // get cell scales
            //~ scalar yOrtho;
            //~ if (useYEff_)
            //~ {
                //~ yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
            //~ }
            //~ else
            //~ {
                //~ yOrtho = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
            //~ }
            //~ scalar V = mesh_.V()[cellI];
            //~ scalar l = Foam::pow(V, 0.333);

            // user-defined scaling
            epsilonIB[bCell] = epsilonIB[bCell]*scaleCoeff_;
            GIB[bCell] = GIB[bCell]*scaleCoeff_;

            //~ // old
            //~ if (body_[cellI] > thrSurf_) // yOrtho < l
            //~ {
                //~ epsilonIB[bCell] = epsilonIB[bCell]*l/(yOrtho + l*0.5); // scalingRight
                //~ GIB[bCell] = GIB[bCell]*l/(yOrtho + l*0.5);
            //~ }

            //~ else // yOrtho > l
            //~ {
                //~ epsilonIB[bCell] = epsilonIB[bCell]/l*(yOrtho + l*0.5); // scalingRightInverse
                //~ GIB[bCell] = GIB[bCell]/l*(yOrtho + l*0.5);
            //~ }
        }
    }

    // assign the values for boundary cells
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // assign
        epsilon[cellI] = epsilonIB[bCell];
        G[cellI] = GIB[bCell];
    }

    // calculate maximum epsilon
    scalar inEpsilon = max(epsilonIB); // internal patch fields for walls should be included as well

    // assign the values in in-solid cells
    forAll(surface, cellI)
    {
        if (surface[cellI] == 1.0)
        {
            if (body_[cellI] >= 0.5)
            {
                epsilon[cellI] = inEpsilon;
                G[cellI] = 0.0;
            }
        }
    }

    // correct the inner boudnary cells
    if (cpDisToInner_)
    {
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get cell labels
            label outCellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;
            label inCellI = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

            // assign
            epsilon[inCellI] = epsilon[outCellI];
            G[inCellI] = G[outCellI];
        }
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::createBaseSurface
(
    volScalarField& surface,
    word surfType,
    scalar boundaryVal
)
{
    if (surfType == "setValue" or surfType == "switched")
    {
        ibInterpolation_->setUpSurface(surface, boundaryVal);
    }

    else if (surfType == "lambdaBased")
    {
        ibInterpolation_->setLambdaBasedSurface(surface, boundaryVal);
    }

    else if (surfType == "onlyInnerValue")
    {
        ibInterpolation_->setOnlyInnerSurface(surface, boundaryVal);
    }

    else
    {
        FatalError << "Surface type " << surfType << " not implemented" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::updateSurface
(
    volScalarField& surface,
    word surfType
)
{
    if (surfType == "switched")
    {
        ibInterpolation_->updateSwitchSurface(surface, ibDirichletBCs_->getYPlusi(), ibDirichletBCs_->getYPlusLam());
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctY
(
    volScalarField& y
)
{
    forAll(boundaryCells_[Pstream::myProcNo()], bCell)
    {
        // get cell label
        label cellI = boundaryCells_[Pstream::myProcNo()][bCell].bCell_;

        // assign distance to boundary
        if (useYEff_)
        {
            y[cellI] = boundaryCells_[Pstream::myProcNo()][bCell].yEff_;
        }
        else
        {
            y[cellI] = boundaryCells_[Pstream::myProcNo()][bCell].yOrtho_;
        }
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::cutFInBoundaryCells
(
    volVectorField& f
)
{
    ibInterpolation_->cutFInBoundaryCells(f);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::cutUInBoundaryCells
(
    volVectorField& U
)
{
    ibInterpolation_->cutUInBoundaryCells(U);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::cutPhiInBoundaryCells
(
    surfaceScalarField& phi
)
{
    ibInterpolation_->cutPhiInBoundaryCells(phi);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::enforceUiInBody
(
    volVectorField& U,
    volVectorField& Ui
)
{
    // loop over cells
    forAll(mesh_.C(), cellI)
    {
        if (body_[cellI] < 1.0)
        {
            continue;
        }

        // check if cellI is an inner boundary cell
        bool toCont = false;
        forAll(boundaryCells_[Pstream::myProcNo()], bCell)
        {
            // get the cell label
            label cellB = boundaryCells_[Pstream::myProcNo()][bCell].iCell_;

            // check
            if (cellB == cellI)
            {
                toCont = true;
                break;
            }
        }

        if (not toCont)
        {
            U[cellI] = Ui[cellI];
        }
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::bound
(
    volScalarField& phi,
    dimensionedScalar& phiMin
)
{
    // loop over cells
    forAll(mesh_.C(), cellI)
    {
        if (body_[cellI] >= 0.5)
        {
            continue;
        }

        phi[cellI] = max(phi[cellI], phiMin.value());
    }
}

// ************************************************************************* //
