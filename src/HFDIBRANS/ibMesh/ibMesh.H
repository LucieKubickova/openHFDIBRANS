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

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
class ibMesh
{
    private:

    // Private data types ------------------------------------------

    // mesh related parameters
    const fvMesh& mesh_;
    //~ scalar VAve_;

    // lambda related things
    const volScalarField& body_;
    //~ scalar thrSurf_;
    //~ scalar intSpan_;

    // options
    //~ word boundarySearch_;
    word stlName_;
    //~ bool excludeWalls_;
    //~ word excludePatch_;
    //~ bool readSurfNorm_;
    //~ bool averageV_;
    //~ bool readL_;
    //~ scalar valueL_;
    bool sdBasedLambda_;
    //~ word surfAreaType_;
    //~ bool yFromCutEdges_; // LK: experimental

    // stl things
    word stlPath_;
    autoPtr<triSurfaceMesh> bodySurfMesh_;
    autoPtr<triSurface> triSurf_;
    autoPtr<triSurfaceSearch> triSurfSearch_;

    // saves

    // references to shared data

    // interpolation info lists

    // dictionaries
    IOdictionary HFDIBDEMDict_;
    //~ IOdictionary fvSchemes_;
    //~ dictionary HFDIBInnerSchemes_;

    // private functions
    // choose interpolation function
    //~ autoPtr<ibScheme> chosenInterpFunc(word name);

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

    //~ void findNeighborInBody(label& cellI, scalar threshold, label& iCell, label& iFace, label& iProc, bool& isBoundary);
    //~ bool isWallCell(label& cellI);
    //~ bool isOnPatch(label& cellI);

    // best instance in direction
    //~ label getFaceInDir(label cellI, vector dir);
    //~ label getEdgeInDir(label faceI, label cellI, vector dir);
    //~ label getVertInDir(label edgeI, label cellI, vector dir);

    // surface creation
    //~ void setUpSurface(volScalarField& surface, scalar boundaryVal);
    //~ void setLambdaBasedSurface(volScalarField& surface, scalar boundaryVal);
    //~ void setOnlyInnerSurface(volScalarField& surface, scalar boundaryVal);
    //~ void updateSwitchSurface(volScalarField& surface, volScalarField& yPlusi, scalar yPlusLam);
    // TODO: blended surface

    // surface correction
    //~ void correctSurfaceByNormal(volSymmTensorField& normSurface, volScalarField& surface, scalar bodyOnLimit);

    // post processing of distance
    //~ scalar getCellSize(label cellI);
    
    // work with stl file
    //~ void getClosestPointAndNormal(const point& startPoint, const vector& span, point& closestPoint, vector& normal);
    //~ scalar createCutCellAndSurface(label cellI, vector& normal, point& surfPoint);
    //~ scalar calculateTriangleArea(point p0, point p1, point p2);
    //~ void createCutCellAndCenter(label cellI, vector& surfNorm, point& surfPoint);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "ibMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
