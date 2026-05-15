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

#ifndef ibMesh_H
#define ibMesh_H

#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "SVD.H"
#include "scalarMatrices.H"
#include "OFstream.H"
#include "Pstream.H"

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "triSurfaceMesh.H"

#include "meshSearch.H"
#include "List.H"
#include "ListOps.H"

#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

#include "triSurfaceMesh.H"
#include "triangleFuncs.H"

#include "ibCutCell.H"
#include "intPoint.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
class ibMesh
{
    private:

    // Private data types ------------------------------------------

    // mesh related parameters
    const fvMesh& mesh_;

    // lambda related things
    const volScalarField& body_;

    // options
    word stlName_;
    bool readL_;
    scalar valueL_;
    bool sdBasedLambda_;
    word cutCellType_;

    // stl things
    word stlPath_;
    autoPtr<triSurfaceMesh> bodySurfMesh_;
    autoPtr<triSurface> triSurf_;
    autoPtr<triSurfaceSearch> triSurfSearch_;

    // dictionaries
    IOdictionary HFDIBDEMDict_;

    public:

    // Constructors ----------------------------------------------------
    ibMesh
    (
        const fvMesh& mesh,
        const volScalarField& body
    );

    // Destructors -----------------------------------------------------
    ~ibMesh();

    // Public methods --------------------------------------------------
    // is point in cell
    bool pointInCell(point pToCheck, label cToCheck);

    // patch checks
    bool isWallCell(label& cellI);
    bool isOnPatch(label& cellI, word& patchName);

    // best instance in direction
    label getFaceInDir(label& cellI, vector& dir);
    label getEdgeInDir(label& faceI, label& cellI, vector& dir);
    label getVertInDir(label& edgeI, label& cellI, vector& dir);

    vector getClosestPoint
    (
        vector ibPoint,
        intPoint& cPoint
    );

    // work with stl files
    void getClosestPointAndNormal(const point& startPoint, const vector& span, point& closestPoint, vector& normal);
    scalar createCutCellAndSurface(label cellI, vector& normal, point& surfPoint);
    scalar calculateTriangleArea(point p0, point p1, point p2);
    void createCutCellAndCenter(label cellI, vector& surfNorm, point& surfPoint);

    // post processing of distance
    scalar getCellSize(label cellI);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
