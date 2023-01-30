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
    DynamicList<label>& boundaryFaces,
    List<bool>& isWallCell
)
:
mesh_(mesh),
boundaryCells_(boundaryCells),
boundaryDists_(boundaryDists),
boundaryFaces_(boundaryFaces),
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
kappa_(0.41),
E_(9.8),
Cmu_(0.09),
Ceps2_(1.9),
beta1_(0.075)
{
    // prepare list to save yPlus in boundary cells
    yPlusi_.setSize(boundaryCells_.size());
    ULogScales_.setSize(boundaryCells_.size());

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
void ibDirichletBCs::pAtIB
(
    List<scalar>& pIB,
    word BCType,
    const volScalarField& p
)
{
    if (BCType == "zeroGradient")
    {
        forAll(boundaryCells_, bCell)
        {
            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // copy pressure from outer cell to the inner cell
            pIB[bCell] = p[cellI];
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
            yPlusi_[bCell] = yPlus;
            kLogScales_[bCell] = uTau/nu[cellI];

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

        forAll(boundaryCells_, bCell)
        {
            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // get distance to the surface
            scalar ds = boundaryDists_[bCell].first();

            // compute magnitude of snGrad of U at the surface
            vector zeroU = vector::zero;
            vector snGradU = (zeroU - U[cellI])/ds;
            scalar magGradUWall = mag(snGradU);

            // compute local Reynolds number
            scalar Rey = ds*Foam::sqrt(k[cellI])/nu[cellI];
            
            // compute normalized variables
            const scalar yPlus = Cmu25_*Rey;
            const scalar uPlus = (1/kappa_)*Foam::log(E_*yPlus);

            // save yPlus
            yPlusi_[bCell] = yPlus;

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

                omegaIB[bCell] = lamFrac*omegaVis + turbFrac*omegaLog;
                GIB[bCell] = lamFrac*G[cellI] + turbFrac*sqr(uStar*magGradUWall*ds/uPlus)/(nu[cellI]*kappa_*yPlus);
            }

            else
            {
                if (yPlus < yPlusLam_)
                {
                    omegaIB[bCell] = 6*nu[cellI]/(beta1_*Foam::sqr(ds));
                    GIB[bCell] = G[cellI];
                }

                else
                {
                    const scalar uStar = Foam::sqrt(Cmu5_*k[cellI]);

                    omegaIB[bCell] = uStar/(Cmu5_*kappa_*ds);
                    GIB[bCell] = sqr(uStar*magGradUWall*ds/uPlus)/(nu[cellI]*kappa_*yPlus);
                }
            }

            if (isWallCell_[bCell])
            {
                omegaIB[bCell] *= 2.0; // provisorial solution
                GIB[bCell] *= 2.0;
            }
        }
    }

    else
    {
        FatalError << omegaWF_ << " condition for omega and G not implemented at the IB" << exit(FatalError);
    }
}

//---------------------------------------------------------------------------//
void ibDirichletBCs::phiAtIB
(
    List<scalar>& phiIB,
    word BCType
)
{
    if (BCType == "noFlux")
    {
        forAll(boundaryFaces_, bFace)
        {
            phiIB[bFace] = 0.0;
        }
    }

    else
    {
        FatalError << BCType << " condition for flux not implemented at the IB" << exit(FatalError);
    }
}

// ************************************************************************* //
