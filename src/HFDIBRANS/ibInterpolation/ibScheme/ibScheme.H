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

#ifndef ibScheme_H
#define ibScheme_H

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
class ibScheme
{
    public:

    // Constructors ----------------------------------------------------
    ibScheme();

    // Destructors -----------------------------------------------------
    virtual ~ibScheme();

    // Public methods --------------------------------------------------
    // -- virtual wrapper functions for interpolation
    virtual scalar interpolate
    (
        volScalarField& phi,
        List<scalar>& phiPs,
        scalar& dirichletVal,
        scalar& scale,
        scalar& ds,
        List<intPoint>& intInfo,
        label& cellI
    ) = 0;

    virtual vector interpolate
    (
        volVectorField& phi,
        List<vector>& phiPs,
        vector& dirichletVal,
        scalar& scale,
        scalar& ds,
        List<intPoint>& intInfo,
        label& cellI
    ) = 0;

    virtual tensor interpolate
    (
        volTensorField& phi,
        List<tensor>& phiPs,
        tensor& dirichletVal,
        scalar& scale,
        scalar& ds,
        List<intPoint>& intInfo,
        label& cellI
    ) = 0;

    // get interpolation order
    label getIntOrder
    (
        List<intPoint>& intPoints
    );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
