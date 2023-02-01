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

#include "ibDirichletBCs.H"

using namespace Foam;

//---------------------------------------------------------------------------//
ibDirichletBCs::ibDirichletBCs
(
    const fvMesh& mesh,
    word simulationType,
    DynamicList<Tuple2<label,label>>& boundaryCells,
    List<Tuple2<scalar,scalar>>& boundaryDists,
    List<Tuple2<bool,label>>& isWallCell
)
:
mesh_(mesh),
boundaryCells_(boundaryCells),
boundaryDists_(boundaryDists),
isWallCell_(isWallCell),
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
simulationType_(simulationType),
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
kappa_(0.41),
E_(9.8),
Cmu_(0.09),
Ceps2_(1.9),
beta1_(0.075)
{
    // ONLY FOR KOMEGA MODELS FOR NOW
    HFDIBBCsDict_ = HFDIBDEMDict_.subDict("wallFunctions");
    HFDIBBCsDict_.lookup("k") >> kWF_;
    HFDIBBCsDict_.lookup("omega") >> omegaWF_;

    // compute turbulence parameters
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
void ibDirichletBCs::UAtIB
(
    List<vector>& UIB,
    word BCType
)
{
    if (simulationType_ == "laminar" or BCType == "noSlip")
    {
        forAll(UIB, bCell)
        {
            UIB[bCell] = ibZero(UIB[bCell]);
        }
    }

    else
    {
        FatalError << BCType << " condition for U in " << simulationType_ << " not implemented at the IB" << exit(FatalError);
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
        forAll(boundaryCells_, bCell)
        {
            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // get distance to the surface
            scalar ds = boundaryDists_[bCell].first();

            // compute the friction velocity
            scalar uTau = Cmu25_*Foam::sqrt(k[cellI]);

            // compute yPlus
            scalar yPlus = uTau*ds/nu[cellI];

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
            kIB[bCell] = max(kIB[bCell], small);
        }
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
        nearWallDist yWall(mesh_);

        forAll(boundaryCells_, bCell)
        {
            // reset fields
            omegaIB[bCell] *= 0.0;
            GIB[bCell] *= 0.0;

            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // prepare list of distances and weights
            List<scalar> distances;
            List<scalar> weights;

            if (isWallCell_[bCell].first())
            {
                // set size of dss
                distances.setSize(2);
                weights.setSize(2);

                // get face label
                label faceI = isWallCell_[bCell].second();

                // get patch label
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                // get local face label
                label lfaceI = mesh_.boundaryMesh()[patchI].whichFace(faceI);

                // get near wall distance
                distances[1] = yWall[patchI][lfaceI];
            }

            else
            {
                // set size of dss
                distances.setSize(1);
                weights.setSize(1);
            }

            // get distance to the surface
            distances[0] = boundaryDists_[bCell].first();

            // calculate weights
            forAll(weights, wI)
            {
                weights[wI] = 1.0/weights.size();
            }

            // loop over walls and surfaces
            forAll(distances, dsI)
            {
                // get distance and weight
                scalar ds = distances[dsI];
                scalar w = weights[dsI];

                // compute magnitude of snGrad of U at the surface
                vector zeroU = vector::zero;
                vector snGradU = (zeroU - U[cellI])/ds;
                scalar magGradUWall = mag(snGradU);

                // compute local Reynolds number
                scalar Rey = ds*Foam::sqrt(k[cellI])/nu[cellI];
                
                // compute normalized variables
                const scalar yPlus = Cmu25_*Rey;
                const scalar uPlus = (1/kappa_)*Foam::log(E_*yPlus);

                // compute the values at the surface
                if (blended)
                {
                    const scalar lamFrac = Foam::exp(-Rey/11);
                    const scalar turbFrac = 1 - lamFrac;

                    const scalar uStar = Foam::sqrt
                    (
                        lamFrac*nu[cellI]*magGradUWall + turbFrac*Cmu5_*k[cellI]
                    );

                    const scalar omegaVis = 6*nu[cellI]/(beta1_*Foam::sqr(ds));
                    const scalar omegaLog = uStar/(Cmu5_*kappa_*ds);

                    omegaIB[bCell] += w*(lamFrac*omegaVis + turbFrac*omegaLog);
                    GIB[bCell] += w*(lamFrac*G[cellI] + turbFrac*sqr(uStar*magGradUWall*ds/uPlus)/(nu[cellI]*kappa_*yPlus));
                }

                else
                {
                    if (yPlus < yPlusLam_)
                    {
                        omegaIB[bCell] += w*(6*nu[cellI]/(beta1_*Foam::sqr(ds)));
                        GIB[bCell] += w*(G[cellI]);
                    }

                    else
                    {
                        const scalar uStar = Foam::sqrt(Cmu5_*k[cellI]);

                        omegaIB[bCell] += w*(uStar/(Cmu5_*kappa_*ds));
                        GIB[bCell] += w*(sqr(uStar*magGradUWall*ds/uPlus)/(nu[cellI]*kappa_*yPlus));
                    }
                }
            }
        }
    }

    else
    {
        FatalError << omegaWF_ << " condition for omega and G not implemented at the IB" << exit(FatalError);
    }
}

// ************************************************************************* //
