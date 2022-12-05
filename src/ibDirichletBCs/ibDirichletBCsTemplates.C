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

//~ #include "ibDirichletBCs.H"

using namespace Foam;

//---------------------------------------------------------------------------//
template<typename Type>
void ibDirichletBCs::UAtIB
(
    List<Type>& UIB,
    word BCType,
    const volScalarField& k,
    const volScalarField& nu
)
{
    if (simulationType_ == "laminar" or BCType == "noSlip")
    {
        forAll(UIB, bCell)
        {
            UIB[bCell] = ibZero(UIB[bCell]);
        }
    }

    else if (simulationType_ == "HFDIBRAS" and BCType == "slip")
    {
        forAll(boundaryCells_, bCell)
        {
            // get cell label
            label cellI = boundaryCells_[bCell].first();

            // get distance to the surface
            scalar ds = boundaryDists_[bCell].first();

            // compute friction velocity
            scalar uTau = Cmu25_*Foam::sqrt(k[cellI]);

            // compute yPlus
            scalar yPlus = uTau*ds/nu[cellI];

            // saves for later interpolation
            yPlusi_[bCell] = yPlus;
            ULogScales_[bCell] = E_*uTau/nu[cellI];

            // compute the value at the surface
            if (yPlus > yPlusLam_)
            {
                UIB[bCell] = 1/kappa_*Foam::log(E_*yPlus)*uTau*ibOne(UIB[bCell]);
            }

            else
            {
                UIB[bCell] = yPlus*uTau*ibOne(UIB[bCell]);
            }
        }
    }

    else
    {
        FatalError << "Dirichlet BCs at the IB for " << simulationType_ << " not implemented" << exit(FatalError);
    }
}

// ************************************************************************* //
