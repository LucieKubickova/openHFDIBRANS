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
),
ibInterpolation_(mesh, body, boundaryCells_, boundaryDists_, isBoundaryCell_),
ibDirichletBCs_(mesh, body, boundaryCells_, boundaryDists_, isBoundaryCell_)
{
    // read HFDIBDEM dictionary
    save_ = HFDIBDEMDict_.lookupOrDefault<bool>("saveIntInfo", false);
    cpOmegaToInner_ = HFDIBDEMDict_.lookupOrDefault<bool>("copyOmegaToInner", false);
    scaleDisG_ = HFDIBDEMDict_.lookupOrDefault<bool>("scaleDisG", false);

    // read fvSchemes
    HFDIBOuterSchemes_ = fvSchemes_.subDict("HFDIBSchemes").subDict("outerSchemes");

    // identify boundary cells
    ibInterpolation_.findBoundaryCells();
    if (save_)
    {
        // save boundary cells as cell sets
        ibInterpolation_.saveBoundaryCells();
    }

    // set size to lists
    ibDirichletBCs_.setSizeToLists();

    // compute distance to the immersed boundary
    boundaryDists_.setSize(boundaryCells_.size());
    ibInterpolation_.calculateDistToBoundary();

    // calculate interpolation points
    ibInterpolation_.calculateInterpolationPoints();

    // save data
    if (save_)
    {
        // create output directory to save data for python
        word outDir = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/ZZ_python";
        if (!isDir(outDir))
        {
            mkDir(outDir);
        }

        // save interpolation data
        ibInterpolation_.saveInterpolationInfo(outDir, "interpolationInfo.dat");
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
    volVectorField& Ui
)
{
    // reset imposed field
    Ui *= 0.0;

    // calculate values at the immersed boundary
    List<vector> UIB;
    UIB.setSize(boundaryCells_.size());
    ibDirichletBCs_.UAtIB(UIB, "noSlip");

    // get references
    volScalarField& yPlusi = ibDirichletBCs_.getYPlusi();
    scalar yPlusLam = ibDirichletBCs_.getYPlusLam();
        
    // calculate log scales for interpolation
    List<scalar> logScales;
    logScales.setSize(boundaryCells_.size());

    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // get distance to surface
        scalar ds = boundaryDists_[bCell].first();

        // calculate the local log scale
        logScales[bCell] = yPlusi[cellI]/ds*ibDirichletBCs_.getE();
    }

    // read interpolation schemes from fvSchemes
    ITstream UIBScheme = HFDIBOuterSchemes_.lookup("U");
    word interpType = UIBScheme[0].wordToken();

    // interpolation
    if (interpType == "unifunctional")
    {
        ibInterpolation_.unifunctionalInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales);
    }

    else if (interpType == "switched")
    {
        ibInterpolation_.switchedInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "outerInner")
    {
        ibInterpolation_.outerInnerInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "inner")
    {
        ibInterpolation_.innerInterp<vector, volVectorField>(UIBScheme, U, Ui, UIB, logScales, yPlusi, yPlusLam);
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
    kIB.setSize(boundaryCells_.size());

    // compute values at the immersed boundary
    ibDirichletBCs_.kAtIB(kIB, k, nu);

    // get references
    volScalarField& yPlusi = ibDirichletBCs_.getYPlusi();
    scalar yPlusLam = ibDirichletBCs_.getYPlusLam();

    // calculate log scales for interpolation
    List<scalar> logScales;
    logScales.setSize(boundaryCells_.size());

    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // get distance to surface
        scalar ds = boundaryDists_[bCell].first();

        // calculate the local log scale
        logScales[bCell] = yPlusi[cellI]/ds;
    }

    // read interpolation schemes from fvSchemes
    ITstream kIBScheme = HFDIBOuterSchemes_.lookup("k");
    word interpType = kIBScheme[0].wordToken();

    // interpolation
    if (interpType == "unifunctional")
    {
        ibInterpolation_.unifunctionalInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales);
    }

    else if (interpType == "switched")
    {
        ibInterpolation_.switchedInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "outerInner")
    {
        ibInterpolation_.outerInnerInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales, yPlusi, yPlusLam);
    }

    else if (interpType == "inner")
    {
        ibInterpolation_.innerInterp<scalar, volScalarField>(kIBScheme, k, ki, kIB, logScales, yPlusi, yPlusLam);
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
void openHFDIBRANS::updateUTau
(
    volScalarField& k
)
{
    ibDirichletBCs_.updateUTauAtIB(k);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctNut
(
    volScalarField& k,
    volScalarField& nu
)
{
    ibDirichletBCs_.correctNutAtIB(k, nu);
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

    omegaIB.setSize(boundaryCells_.size());
    GIB.setSize(boundaryCells_.size());

    // calculate values at the immersed boundary
    ibDirichletBCs_.omegaGAtIB(omegaIB, GIB, G, U, k, nu);

    // omega scaling
    if (scaleDisG_)
    {
        forAll(boundaryCells_, bCell)
        {
            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // get cell scales
            scalar yOrtho = boundaryDists_[bCell].first();
            scalar V = mesh_.V()[cellI];
            scalar l = Foam::pow(V, 0.333);

            // assign
            //~ omegaIB[bCell] = omegaIB[bCell]/yOrtho*l; // scaling
            //~ GIB[bCell] = GIB[bCell]/yOrtho*l;
            omegaIB[bCell] = omegaIB[bCell]*l/(yOrtho + l*0.5); // scalingRight
            GIB[bCell] = GIB[bCell]*l/(yOrtho + l*0.5);
            //~ omegaIB[bCell] = omegaIB[bCell]/l*(yOrtho + l*0.5); // scalingRightInverse
            //~ GIB[bCell] = GIB[bCell]/l*(yOrtho + l*0.5);
        }
    }

    // assign the values for boundary cells
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

        // assign
        omega[cellI] = omegaIB[bCell];
        G[cellI] = GIB[bCell];
    }

    // calculate maximum omega
    scalar inOmega = max(omegaIB); // internal patch fields for walls should be included as well

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

    // correct the inner boundary cells
    if (cpOmegaToInner_)
    {
        forAll(boundaryCells_, bCell)
        {
            // get cell labels
            label outCellI = boundaryCells_[bCell].first();
            label inCellI = boundaryCells_[bCell].second();

            // assign
            omega[inCellI] = omega[outCellI];
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

    epsilonIB.setSize(boundaryCells_.size());
    GIB.setSize(boundaryCells_.size());

    // calculate values at the immersed boundary
    ibDirichletBCs_.epsilonGAtIB(epsilonIB, GIB, G, U, k, nu);

    // epsilon scaling
    if (scaleDisG_)
    {
        forAll(boundaryCells_, bCell)
        {
            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // get cell scales
            scalar yOrtho = boundaryDists_[bCell].first();
            scalar V = mesh_.V()[cellI];
            scalar l = Foam::pow(V, 0.333);

            // assign
            //~ epsilonIB[bCell] = epsilonIB[bCell]/yOrtho*l; // scaling
            //~ GIB[bCell] = GIB[bCell]/yOrtho*l;
            epsilonIB[bCell] = epsilonIB[bCell]*l/(yOrtho + l*0.5); // scalingRight
            GIB[bCell] = GIB[bCell]*l/(yOrtho + l*0.5);
            //~ epsilonIB[bCell] = epsilonIB[bCell]/l*(yOrtho + l*0.5); // scalingRightInverse
            //~ GIB[bCell] = GIB[bCell]/l*(yOrtho + l*0.5);
        }
    }

    // assign the values for boundary cells
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].first();

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
                G[cellI] = small;
            }
        }
    }

    //~ // correct the inner boudnary cells
    //~ forAll(boundaryCells_, bCell)
    //~ {
        //~ // get cell labels
        //~ label outCellI = boundaryCells_[bCell].first();
        //~ label inCellI = boundaryCells_[bCell].second();

        //~ // assign
        //~ epsilon[inCellI] = epsilon[outCellI];
    //~ }
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
        ibInterpolation_.setUpSurface(surface, boundaryVal);
    }

    else if (surfType == "lambdaBased")
    {
        ibInterpolation_.setLambdaBasedSurface(surface, boundaryVal);
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
        ibInterpolation_.updateSwitchSurface(surface, ibDirichletBCs_.getYPlusi(), ibDirichletBCs_.getYPlusLam());
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::cutFInBoundaryCells
(
    volVectorField& f
)
{
    ibInterpolation_.cutFInBoundaryCells(f);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::cutUInBoundaryCells
(
    volVectorField& U
)
{
    ibInterpolation_.cutUInBoundaryCells(U);
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::cutPhiInBoundaryCells
(
    surfaceScalarField& phi
)
{
    ibInterpolation_.cutPhiInBoundaryCells(phi);
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
        forAll(boundaryCells_, bCell)
        {
            // get the cell label
            label cellB = boundaryCells_[bCell].second();

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
