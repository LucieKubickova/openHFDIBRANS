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
    class for interpolation

SourceFiles

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef lineIntInfo_H
#define lineIntInfo_H

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

#include "intPoint.H"

using namespace Foam;

//---------------------------------------------------------------------------//
class lineIntInfo
{
    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

    // mesh reference
    const Foam::fvMesh& mesh_;

    // lists
    List<label> ibCells_;
    List<point> ibPoints_;
    List<vector> ibNormals_;

    private:

    autoPtr<List<List<intPoint>>> intPointsPtr_; // Points used for interp.

    intPoint  findIntPoint
    (
        intPoint& fromP,
        point& endP
    );

    void correctIntPoint
    (
        point ibPoint,
        intPoint& cPoint
    );

    vector getClosestPoint
    (
        vector ibPoint,
        intPoint& cPoint
    );

    label getFaceInDir
    (
        const intPoint& retPoint,
        const label prevFace
    );

    bool pointInCell
    (
        point pToCheck,
        label cToCheck
    );

    void syncIntPoints();

    public:

    // Constructors ------------------------------------------------
    lineIntInfo
    (
        const fvMesh& mesh,
        List<label> ibCells,
        List<point> ibPoints,
        List<vector> ibNormals
    );

    // Destructors -------------------------------------------------
    ~lineIntInfo();

    void setIntpInfo();

    void resetIntpInfo(label newSize)
    {
        intPointsPtr_.reset(new List<List<intPoint>>(newSize));
    }

    List<List<intPoint>>& getIntPoints()
    {
        return intPointsPtr_();
    }

    const List<List<intPoint>>& getIntPoints() const
    {
        return intPointsPtr_();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
