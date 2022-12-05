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

Description
    implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with (level-set) DEM (Isoz and Sourek, 2020)

SourceFiles
    openHFDIBDEM.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef ibInterpolation_H
#define ibInterpolation_H

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

#include "interpolationInfo.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
class ibInterpolation
{
    private:

    // Private data types ------------------------------------------

    // mesh related parameters
    const fvMesh& mesh_;
    scalar VAve_;

    // lambda related things
    const volScalarField& body_;
    volVectorField surfNorm_;
    volVectorField surfTan_;
    scalar thrSurf_;
    scalar intSpan_;

    // references to shared data
    DynamicList<Tuple2<label,label>>& boundaryCells_;
    List<Tuple2<scalar,scalar>>& boundaryDists_;

    // interpolation info lists
    List<interpolationInfo> intInfoList_;

    // dictionaries
    IOdictionary HFDIBDEMDict_;
    dictionary HFDIBInterpDict_;

    // private functions
    // peace-wise interpolation
    template <typename Type, typename volTypeField>
    Type polynom(volTypeField& phi, interpolation<Type>& interpPhi, Type& dirichletVal, label bCell);

    template <typename Type, typename volTypeField>
    Type logarithm(volTypeField& phi, interpolation<Type>& interpPhi, scalar& dirichletVal, scalar& logScale, label bCell);

    public:

    // Constructors ----------------------------------------------------
    ibInterpolation(const fvMesh& mesh, const volScalarField& body, DynamicList<Tuple2<label,label>>& boundaryCells, List<Tuple2<scalar,scalar>>& boundaryDists);

    // Destructors -----------------------------------------------------
    ~ibInterpolation();

    // Public methods --------------------------------------------------
    
    // compute imposed field
    template <typename Type, typename volTypeField>
    void polynomialInterp(volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals);

    template <typename Type, typename volTypeField>
    void logarithmicInterp(volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& logScales);

    template <typename Type, typename volTypeField>
    void polySwitchLogInterp(volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& logScales, List<scalar>& yPlusi, scalar yPlusLam);

    // TODO: blended interpolation
    
    // cell identification
    void findBoundaryCells();
    void createSurface(volScalarField& surface);

    // surface distance calculation
    void calculateDistToBoundary();
    
    // calculate interpolation points
    void calculateInterpolationPoints();

    // calculate surface directions
    void calculateSurfTans(volVectorField& U);
    void calculateSurfNorm();
    void correctSurfNorm();

    // find cell
    Tuple2<vector,Tuple2<label,label>> findCellCustom(vector& prevPoint, label& startCell, label& startProc, vector& gradToBody, scalar& intDist);

    // save
    void saveInterpolationInfo(word outDir, word fileName);
    void saveBoundaryCells();
    void saveCellSet(List<label>& listToSave, word fileName);

    // return the surface threshold
    scalar getSurfaceThreshold()
    {
        return thrSurf_;
    }

    // return the surface normals field
    volVectorField& getSurfNorm()
    {
        return surfNorm_;
    }

    // return the surface tangents field
    volVectorField& getSurfTan()
    {
        return surfTan_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "ibInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
