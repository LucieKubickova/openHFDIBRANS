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

Class
    Foam::ibCutCell

Description
    cutCell class to be used for HFDIBRANS turbulence models

SourceFiles
    ibCutCell.C

Contributors
	Federico Municchi (2016),
	Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),

\*---------------------------------------------------------------------------*/

#ifndef ibCutCell_H
#define ibCutCell_H

#include "fvMesh.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
						Class ibCutCell Declaration
\*---------------------------------------------------------------------------*/

class ibCutCell
{

private:

    // Private Data

		//- fvMesh
		const fvMesh& mesh_;

		//- Surface normal
		const vector& surfNorm_;
		const point& surfPoint_;
		const cell& inCell_;

		DynamicList<DynamicList<point>> v_;

		DynamicList<DynamicList<point>> faces_;
		List<vector> faceCentres_;
		List<vector> faceAreas_;

		//- Cut cell center
		vector C_;

		//- Cut cell volume
		scalar V_;

		//- Input cell center
		vector CIn_;

		//- Input cell volume
		scalar VIn_;

		//- Cut cell type
		// -1 ... in
		//  0 ... cut
		//  1 ... out
		label type_;

		//- Number of cut faces
		label nCutFaces_;

		// Private Member Functions

			//- Make faces
			void makeFaces(const List<point>& vertices);

			//- Make cell center and volume
			void makeCellCentreAndVol();

			//- Make face centers and areas
			void makeFaceCentresAndAreas();

			//- Synchronize normals
			void syncNormals();


public:

    // Constructors

		//- Construct from components
		ibCutCell
		(
			const fvMesh& mesh,
			const vector& surfNorm,
			const point& surfPoint,
			const cell& inCell
		);

    //- Destructor
    ~ibCutCell();

    // Member Functions

		// Access

			//- Get cut cell center
			vector C()
			{
				return C_;
			};

			//- Get cut cell volume
			scalar V()
			{
				return V_;
			};

			//- Get cut cell surface point
			vector S()
			{
				return surfPoint_;
			};

			//- Get input cell center
			vector CIn()
			{
				return CIn_;
			};

			//- Get input cell volume
			scalar VIn()
			{
				return VIn_;
			};

			//- Get cut cell face areas
			List<vector> Sf()
			{
				return faceAreas_;
			};

			//- Get cutting surface normal
			vector Ss()
			{
				return surfNorm_;
			};

			//- Get cut cell faces
			DynamicList<DynamicList<point>> faces()
			{
				return faces_;
			};

			//- Get cut cell type
			label status()
			{
				return type_;
			};

		//- Get orthogonal distance from cut cell center to the wall
		scalar yOrtho();

		//- Get 1/2 of the distance from surface point to the opposite face
		scalar yOrthoEst();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
