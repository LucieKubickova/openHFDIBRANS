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
surface_("surface", body),
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
ibInterpolation_(mesh, body, boundaryCells_, boundaryDists_),
ibDirichletBCs_(mesh, simulationType, boundaryCells_, boundaryDists_)
{
    // read HFDIBDEM dictionary
    HFDIBDEMDict_.lookup("saveIntInfo") >> save_;

    // identify boundary cells
    ibInterpolation_.findBoundaryCells();

    // correct ibDirichletBCs lists
    ibDirichletBCs_.correctLists();

    // create surface field
    ibInterpolation_.createSurface(surface_);

    // compute distance to the surface
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

    else if (simulationType_ == "HFDIBRAS")
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

        // polynomial interpolation for normal velocity component
        ibInterpolation_.polynomialInterp<scalar, volScalarField>(UNorm, UNormi, UNormIB);
        // TODO: blendedInterp

        // construct the Ui field from UTani and UNormi
        Ui = UTani*ibInterpolation_.getSurfTan() + UNormi*ibInterpolation_.getSurfNorm();
    }

    else
    {
        FatalError << "Interpolation for " << simulationType_ << " not implemented" << exit(FatalError);
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
    List<scalar>& logScales = ibDirichletBCs_.getkLogScales();
    List<scalar>& yPlusi = ibDirichletBCs_.getYPlusi();
    scalar yPlusLam = ibDirichletBCs_.getYPlusLam();

    // switch interpolation
    ibInterpolation_.polySwitchLogInterp<scalar, volScalarField>(k, ki, kIB, logScales, yPlusi, yPlusLam);
    // TODO: blendedInterp
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

    // assign the values in in-solid cells
    scalar inOmega = max(omega).value();

    forAll(surface_, cellI)
    {
        if (surface_[cellI] == 1.0)
        {
            if (body_[cellI] >= 0.5)
            {
                omega[cellI] = inOmega;
                G[cellI] = 0.0;
            }
        }
    }
}

// ************************************************************************* //
