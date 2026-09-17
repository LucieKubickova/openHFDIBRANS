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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Get zero scalar
inline scalar ibDirichletBCs::ibZero(scalar dummy)
{
	return 0.0;
}

//- Get zero vector
inline vector ibDirichletBCs::ibZero(vector dummy)
{
	return vector::zero;
}

//- Get unit scalar
inline scalar ibDirichletBCs::ibOne(scalar dummy)
{
	return 1.0;
}

//- Get unit vector
inline vector ibDirichletBCs::ibOne(vector dummy)
{
	return vector::one;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Get immersed yPlus field
inline volScalarField& ibDirichletBCs::getYPlusi()
{
	return yPlusi_;
}

//- Get yPlusLam
inline scalar ibDirichletBCs::getYPlusLam()
{
	return yPlusLam_;
}

//- Get E coefficient
inline scalar ibDirichletBCs::getE()
{
	return E_;
}

//- Get nut at immersed boundary
inline List<List<scalar>>& ibDirichletBCs::getNutAtIB()
{
	return nutAtIB_;
}

//- Get k at immersed boundary
inline List<List<scalar>>& ibDirichletBCs::getKAtIB()
{
	return kAtIB_;
}

//- Get imposed nut field
inline volScalarField& ibDirichletBCs::getNuti()
{
	return nuti_;
}

//- Get imposed friction velocity field
inline volScalarField& ibDirichletBCs::getUTaui()
{
	return uTaui_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
