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
    List<Tuple2<scalar,scalar>>& boundaryDists
)
:
mesh_(mesh),
boundaryCells_(boundaryCells),
boundaryDists_(boundaryDists),
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
    HFDIBBCsDict_.lookup("nut") >> nutWF_;
    HFDIBBCsDict_.lookup("k") >> kWF_;
    HFDIBBCsDict_.lookup("omega") >> omegaWF_;
    HFDIBBCsDict_.lookup("epsilon") >> epsilonWF_;

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
void ibDirichletBCs::setSizeToLists
(
)
{
    // set size to nutAtIB
    nutAtIB_.setSize(boundaryCells_.size());
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
void ibDirichletBCs::correctNutAtIB
(
    volScalarField& k,
    volScalarField& nu
)
{
    if (nutWF_ == "nutkWallFunction")
    {
        // loop over boundary cells
        forAll(boundaryCells_, bCell)
        {
            // reset field
            nutAtIB_[bCell] *= 0.0;

            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // get distance to the surface
            scalar ds = boundaryDists_[bCell].first();

            // compute the friction velocity
            scalar uTau = Cmu25_*Foam::sqrt(k[cellI]);

            // compute yPlus
            scalar yPlus = uTau*ds/nu[cellI];

            // compute the values at the surface
            if (yPlus > yPlusLam_)
            {
                nutAtIB_[bCell] = nu[cellI]*(yPlus*kappa_/Foam::log(E_*yPlus) - 1.0);
            }
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
        nearWallDist yWall(mesh_);

        forAll(boundaryCells_, bCell)
        {
            // reset fields
            epsilonIB[bCell] *= 0.0;
            GIB[bCell] *= 0.0;

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

            if (yPlus > yPlusLam_)
            {
                epsilonIB[bCell] = Cmu75_*Foam::pow(k[cellI], 1.5)/(kappa_*ds);
                GIB[bCell] = nutAtIB_[bCell] + nu[cellI]*magGradUWall*Cmu25_*Foam::sqrt(k[cellI])/(kappa_*ds);
            }

            else
            {
                epsilonIB[bCell] = 2.0*k[cellI]*nu[cellI]/sqr(ds);
                GIB[bCell] = G[cellI]; // NOTE: not sure about this
            }
        }
    }

    else
    {
        FatalError << epsilonWF_ << " condition for epsilon and G not implemented at the IB" << exit(FatalError);
    }
}


// ************************************************************************* //
