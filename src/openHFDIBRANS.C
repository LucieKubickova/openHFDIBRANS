/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

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
    const volScalarField& body,
    word simulationType
) :
mesh_(mesh),
body_(body),
outSurface_("outSurface", body),
inSurface_("inSurface", body),
simulationType_(simulationType),
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
ibInterpolation_(mesh, body, boundaryCells_, boundaryDists_, boundaryFaces_, isWallCell_),
ibDirichletBCs_(mesh, simulationType, boundaryCells_, boundaryDists_, boundaryFaces_, isWallCell_)
{
    // read HFDIBDEM dictionary
    HFDIBDEMDict_.lookup("saveIntInfo") >> save_;

    // interpolation types
    HFDIBDEMDict_.lookup("UInterpType") >> UInterpType_;

    // simple algorithm settings
    dictionary HFDIBSimpleDict = HFDIBDEMDict_.subDict("simple");
    tolUEqn_ = readScalar(HFDIBSimpleDict.lookup("tolUEqn"));
    tolKEqn_ = readScalar(HFDIBSimpleDict.lookup("tolKEqn"));
    maxUEqnIters_ = readScalar(HFDIBSimpleDict.lookup("maxUEqnIters"));
    maxPEqnIters_ = readScalar(HFDIBSimpleDict.lookup("maxPEqnIters"));
    maxKEqnIters_ = readScalar(HFDIBSimpleDict.lookup("maxKEqnIters"));

    // identify boundary cells
    ibInterpolation_.findBoundaryCells();

    // identify boundary faces
    ibInterpolation_.findBoundaryFaces();

    // check whether boundary cells are adjecent to regular walls
    isWallCell_.setSize(boundaryCells_.size());
    ibInterpolation_.areWallCells();

    // correct ibDirichletBCs lists
    ibDirichletBCs_.correctLists();

    // create surface fields
    ibInterpolation_.createOuterSurface(outSurface_);
    ibInterpolation_.createInnerSurface(inSurface_);

    // compute distance to the immersed boundary
    boundaryDists_.setSize(boundaryCells_.size());
    ibInterpolation_.calculateDistToBoundary();

    // calculate interpolation points
    ibInterpolation_.calculateInterpolationPoints();

    // correct the surface normal
    ibInterpolation_.correctSurfNorm();

    // save data
    if (save_)
    {
        // save boundary cells as cell sets
        ibInterpolation_.saveBoundaryCells();

        // create output directory to save data
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
    volVectorField& Ui,
    const volScalarField& k,
    const volScalarField& nu
)
{
    // reset imposed field
    Ui *= 0.0;

    // interpolate to boundary cells centers
    if (simulationType_ == "laminar")
    {
        // calculate values at the immersed boundary
        List<vector> UIB;
        UIB.setSize(boundaryCells_.size());
        ibDirichletBCs_.UAtIB<vector>(UIB, "noSlip", k, nu);

        // interpolate
        ibInterpolation_.polynomialInterp<vector, volVectorField>(U, Ui, UIB);
    }

    else if (simulationType_ == "HFDIBRAS" and UInterpType_ == "outer")
    {
        // calculate surface tangents
        ibInterpolation_.calculateSurfTans(U);

        // create normal and tangential velocity fields
        volScalarField UTan("UTan", U & ibInterpolation_.getSurfTan());
        volScalarField UNorm("UNorm", U & ibInterpolation_.getSurfNorm());

        // create immposed variants
        volScalarField UTani(UTan*0.0);
        volScalarField UNormi(UNorm*0.0);

        // calculate values at the immersed boundary
        List<scalar> UTanIB;
        List<scalar> UNormIB;

        UTanIB.setSize(boundaryCells_.size());
        UNormIB.setSize(boundaryCells_.size());

        ibDirichletBCs_.UAtIB<scalar>(UTanIB, "slip", k, nu);
        ibDirichletBCs_.UAtIB<scalar>(UNormIB, "noSlip", k, nu);

        // get references
        List<scalar>& logScales = ibDirichletBCs_.getULogScales();
        List<scalar>& yPlusi = ibDirichletBCs_.getYPlusi();
        scalar yPlusLam = ibDirichletBCs_.getYPlusLam();

        // switch interpolation only for tangential velocity component
        ibInterpolation_.polySwitchLogInterp<scalar, volScalarField>(UTan, UTani, UTanIB, logScales, yPlusi, yPlusLam);
        //~ ibInterpolation_.polySwitchConsInterp<scalar, volScalarField>(UTan, UTani, UTanIB, logScales, yPlusi, yPlusLam); // nope, does not work
        //~ ibInterpolation_.noInterp<scalar, volScalarField>(UTan, UTani, UTanIB); // nope, blows up

        // polynomial interpolation for normal velocity component
        ibInterpolation_.polynomialInterp<scalar, volScalarField>(UNorm, UNormi, UNormIB);

        // TODO: blended interpolation
        // TODO: interpolation inside?

        // construct the Ui field from UTani and UNormi
        Ui = UTani*ibInterpolation_.getSurfTan() + UNormi*ibInterpolation_.getSurfNorm();

        // correct boundary cells adjecent to regular walls
        forAll(boundaryCells_, bCell)
        {
            // get cell label
            label cellI = boundaryCells_[bCell].first();

            if (isWallCell_[bCell])
            {
                Ui[cellI] = vector::zero; // provisorial solution
            }
        }
    }

    else if (simulationType_ == "HFDIBRAS" and UInterpType_ == "inner")
    {
        // calculate surface tangents
        ibInterpolation_.calculateSurfTans(U);

        // calculate values at the immersed boundary
        List<scalar> UTanIB;
        List<scalar> UNormIB;

        UTanIB.setSize(boundaryCells_.size());
        UNormIB.setSize(boundaryCells_.size());

        ibDirichletBCs_.UAtIB<scalar>(UTanIB, "slip", k, nu);
        ibDirichletBCs_.UAtIB<scalar>(UNormIB, "noSlip", k, nu);

        // construct UIB from normal and tangential components
        List<vector> UIB;
        UIB.setSize(boundaryCells_.size());

        forAll(boundaryCells_, bCell)
        {
            // get the cell label
            label cellI = boundaryCells_[bCell].first();

            // get normal and tangential direction vectors
            vector tan = ibInterpolation_.getSurfTan()[cellI];
            vector norm = ibInterpolation_.getSurfNorm()[cellI];

            // construct the value at the boundary
            //~ UIB[bCell] = UTanIB[bCell]*tan + UNormIB[bCell]*norm;
            UIB[bCell] = vector::zero;
        }

        // copy values of U to Ui in the outer cells
        ibInterpolation_.noInterp<vector, volVectorField>(U, Ui, UIB);

        // linear interpolation to inner cells
        ibInterpolation_.linearInInterp<vector, volVectorField>(U, Ui, UIB);
    }

    else
    {
        FatalError << "Interpolation for " << simulationType_ << " not implemented" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctP
(
    volScalarField& p
)
{
    // prepare lists
    List<scalar> pIB;
    pIB.setSize(boundaryCells_.size());

    // compute values at the immersed boundary
    ibDirichletBCs_.pAtIB(pIB, "zeroGradient", p);

    // correct pressure
    forAll(boundaryCells_, bCell)
    {
        // get cell label
        label cellI = boundaryCells_[bCell].second();

        // assign pressure value
        p[cellI] = pIB[bCell];
    }
}

//---------------------------------------------------------------------------//
scalar openHFDIBRANS::maxErrorInPBC
(
    volScalarField& p
)
{
    // prepare Lists
    List<scalar> errorInP;
    errorInP.setSize(boundaryCells_.size());

    // compute error in pressure boundary condition
    forAll(boundaryCells_, bCell)
    {
        // get label of the outer cell
        label outerCellI = boundaryCells_[bCell].first();

        // get label of the inner cell
        label innerCellI = boundaryCells_[bCell].second();

        // evaluate the error
        errorInP[bCell] = mag(p[outerCellI] - p[innerCellI]);
    }

    // evaluate the max value
    scalar maxError = max(errorInP);

    // return
    return maxError;
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
    List<scalar>& logScales = ibDirichletBCs_.getkLogScales();
    List<scalar>& yPlusi = ibDirichletBCs_.getYPlusi();
    scalar yPlusLam = ibDirichletBCs_.getYPlusLam();

    // switch interpolation
    ibInterpolation_.polySwitchLogInterp<scalar, volScalarField>(k, ki, kIB, logScales, yPlusi, yPlusLam);
    //~ ibInterpolation_.polySwitchConsInterp<scalar, volScalarField>(k, ki, kIB, logScales, yPlusi, yPlusLam); // nope, does not work

    // TODO: blended interpolation
    // TODO: interpolate inside?
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctOmegaG
(
    volScalarField& omega,
    volScalarField::Internal& G,
    const volVectorField& U,
    volScalarField& k,
    volScalarField& nu
)
{
    // prepare lists
    List<scalar> omegaIB;
    List<scalar> GIB;

    omegaIB.setSize(boundaryCells_.size());
    GIB.setSize(boundaryCells_.size());

    // calculate values at the immersed boundary
    ibDirichletBCs_.omegaGAtIB(omegaIB, GIB, G, U, k, nu);

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
    forAll(outSurface_, cellI)
    {
        if (outSurface_[cellI] == 1.0)
        {
            if (body_[cellI] >= 0.5)
            {
                omega[cellI] = inOmega;
                G[cellI] = 0.0;
            }
        }
    }
}

//---------------------------------------------------------------------------//
void openHFDIBRANS::correctPhi
(
    surfaceScalarField& phi
)
{
    // prepare lists
    List<scalar> phiIB;
    phiIB.setSize(boundaryFaces_.size());

    // calculate values at the immersed boundary
    ibDirichletBCs_.phiAtIB(phiIB, "noFlux");

    // assign values to faces forming the immersed boundary
    forAll(boundaryFaces_, bFace)
    {
        // get face label
        label faceI = boundaryFaces_[bFace];

        // assign the values
        phi[faceI] = phiIB[bFace];
    }
}

// ************************************************************************* //
