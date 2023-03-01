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

#include "ibScheme.H"
#include "constantScheme.H"
#include "linearScheme.H"
#include "quadraticScheme.H"
#include "logarithmicScheme.H"

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
    scalar thrSurf_;
    scalar intSpan_;

    // references to shared data
    DynamicList<Tuple2<label,label>>& boundaryCells_;
    List<Tuple2<scalar,scalar>>& boundaryDists_;
    List<Tuple2<bool,label>>& isWallCell_;

    // interpolation info lists
    List<interpolationInfo> intInfoList_;

    // dictionaries
    IOdictionary HFDIBDEMDict_;
    dictionary HFDIBInterpDict_;

    // private functions
    // choose interpolation function
    autoPtr<ibScheme> chosenInterpFunc(word name);

    public:

    // Constructors ----------------------------------------------------
    ibInterpolation(const fvMesh& mesh, const volScalarField& body, DynamicList<Tuple2<label,label>>& boundaryCells, List<Tuple2<scalar,scalar>>& boundaryDists, List<Tuple2<bool,label>>& isWallCell);

    // Destructors -----------------------------------------------------
    ~ibInterpolation();

    // Public methods --------------------------------------------------
    // -- interpolation approaches
    // unifunctional interpolation
    template <typename Type, typename volTypeField>
    void unifunctionalInterp(ITstream& interpScheme, volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& scales);

    // switched, decision based on local yPlus value
    template <typename Type, typename volTypeField>
    void switchedInterp(ITstream& interpScheme, volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& scales, volScalarField& yPlusi, scalar yPlusLam);
    
    // TODO: blended interpolation
    
    // identification
    void findBoundaryCells();
    void areWallCells();

    // surface creation
    void setUpSurface(volScalarField& surface, scalar boundaryVal);
    void updateSwitchSurface(volScalarField& surface, volScalarField& yPlusi, scalar yPlusLam);
    // TODO: blended surface

    // surface distance calculation
    void calculateDistToBoundary();
    
    // calculate interpolation points
    void calculateInterpolationPoints();

    // surface normals
    void calculateSurfNorm();
    void correctSurfNorm();

    // find cell
    Tuple2<vector,Tuple2<label,label>> findCellCustom(vector& prevPoint, label& startCell, label& startProc, vector& gradToBody, scalar& intDist);

    // save
    void saveInterpolationInfo(word outDir, word fileName);
    void saveBoundaryCells();
    void saveCellSet(List<label>& listToSave, word fileName);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "ibInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
