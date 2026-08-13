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

Application
    generateLambda

Description
	A utility application to generate the lambda field based
	on a provided STL file

\*---------------------------------------------------------------------------*/

#include "IOdictionary.H"
#include "IOobject.H"
#include "dimensionedScalarFwd.H"
#include "fvCFD.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "triSurfaceSearch.H"
#include "volFieldsFwd.H"

#include "stlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
	#include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info << "Preparing lambda\n" << endl;

	// Prepare lambda field
	volScalarField lambda
	(
		IOobject
		(
			"lambda",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimless, 0.0)
	);

	Info << "Reading HFDIBDEMDict\n" << endl;

	// Load dictionary
	IOdictionary HFDIBDEMDict
	(
		IOobject
		(
			"HFDIBDEMDict",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	// Read options
	word stlName;
	HFDIBDEMDict.lookup("stlName") >> stlName;
	scalar thrSurf =
		HFDIBDEMDict.lookupOrDefault<scalar>("surfaceThreshold", 1.0);
	scalar intSpan =
		HFDIBDEMDict.lookupOrDefault<scalar>("interfaceSpan", 1.0);

	// Load the STL
	Info << "Reading the " << stlName << " file\n" << endl;
	autoPtr<triSurfaceMesh> surfMesh
	(
		new triSurfaceMesh(
			IOobject
			(
				stlName + ".stl",
				runTime.constant(),
				"triSurface",
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		)
	);

	autoPtr<triSurface> triSurf(new triSurface(surfMesh()));
	autoPtr<triSurfaceSearch> triSurfSearch(new triSurfaceSearch(triSurf()));

	Info << "Generating lambda based on " << stlName << nl << endl;

	// Initialize STL model and generate lambda
	stlModel model(mesh, thrSurf, intSpan, surfMesh, triSurfSearch);
	model.generateLambda(lambda);

	Info << "Writing lambda\n" << endl;

	// Write lambda field
	lambda.write();

	Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
