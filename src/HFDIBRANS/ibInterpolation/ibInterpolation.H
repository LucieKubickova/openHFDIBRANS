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

#ifndef ibInterpolation_H
#define ibInterpolation_H

#include "ibScheme.H"
#include "constantScheme.H"
#include "linearScheme.H"
#include "quadraticScheme.H"
#include "logarithmicScheme.H"
#include "fixedGradientScheme.H"
#include "zeroGradientScheme.H"

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

    // options
    word boundarySearch_;
    bool excludeWalls_;
    bool readSurfNorm_;
    bool aveYOrtho_;
    bool totalYOrthoAve_;
    scalar aveCoeff_;
    label nAveYOrtho_;
    bool averageV_;

    // saves
    volScalarField yOrthoi_;

    // references to shared data
    List<DynamicList<Tuple2<label,label>>>& boundaryCells_;
    List<List<Tuple2<scalar,scalar>>>& boundaryDists_;
    List<DynamicList<label>>& surfaceCells_;
    List<List<scalar>>& surfaceDists_;
    List<DynamicList<label>>& internalCells_;
    labelField& isBoundaryCell_;

    // interpolation info lists
    List<List<interpolationInfo>> intInfoListBoundary_;
    List<List<interpolationInfo>> intInfoListSurface_;

    // dictionaries
    IOdictionary HFDIBDEMDict_;
    IOdictionary fvSchemes_;
    dictionary HFDIBInnerSchemes_;

    // private functions
    // choose interpolation function
    autoPtr<ibScheme> chosenInterpFunc(word name);

    public:

    // Constructors ----------------------------------------------------
    ibInterpolation
    (
        const fvMesh& mesh,
        const volScalarField& body,
        List<DynamicList<Tuple2<label,label>>>& boundaryCells,
        List<List<Tuple2<scalar,scalar>>>& boundaryDists,
        List<DynamicList<label>>& surfaceCells,
        List<List<scalar>>& surfaceDists,
        List<DynamicList<label>>& internalCells,
        labelField& isBoundaryCell
    );

    // Destructors -----------------------------------------------------
    ~ibInterpolation();

    // Public methods --------------------------------------------------
    // -- interpolation approaches
    // unifunctional interpolation
    template <typename Type, typename volTypeField>
    void unifunctionalInterp(ITstream& ibSchemeName, volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& scales);

    // lambda based interpolation (unifuncional)
    template <typename Type, typename volTypeField>
    void lambdaBasedInterp(ITstream& ibSchemeName, volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& scales);

    // switched, decision based on local yPlus value
    template <typename Type, typename volTypeField>
    void switchedInterp(ITstream& ibSchemeName, volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& scales, volScalarField& yPlusi, scalar yPlusLam);
    
    // outer or inner interpolation, decision based on local yPlus value
    template <typename Type, typename volTypeField>
    void outerInnerInterp(ITstream& ibSchemeName, volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& scales, volScalarField& yPlusi, scalar yPlusLam);

    // inner unifunctional interpolation
    template <typename Type, typename volTypeField>
    void innerInterp(ITstream& ibSchemeName, volTypeField& phi, volTypeField& phii, List<Type>& dirichletVals, List<scalar>& scales, volScalarField& yPlusi, scalar yPlusLam);

    // TODO: blended interpolation
    
    // identification
    void findBoundaryCells();
    void findSurfaceCells();

    // surface creation
    void setUpSurface(volScalarField& surface, scalar boundaryVal);
    void setLambdaBasedSurface(volScalarField& surface, scalar boundaryVal);
    void updateSwitchSurface(volScalarField& surface, volScalarField& yPlusi, scalar yPlusLam);
    // TODO: blended surface

    // surface distance calculation
    void calculateBoundaryDist();
    void calculateSurfaceDist();
    
    // calculate interpolation points
    void calculateInterpolationPoints();
    void getInterpolationPoint(label cellI, point surfPoint, vector surfNormToSend, interpolationInfo& intInfo);
    void setInterpolationOrder(List<interpolationInfo>& intInfoList);

    // surface normals
    void calculateSurfNorm();

    // find cell
    Tuple2<vector,Tuple2<label,label>> findCellCustom(vector& prevPoint, label& startCell, label& startProc, vector& gradToBody, scalar& intDist);

    // save
    void saveInterpolationInfo(word outDir, word fileName);
    void saveBoundaryCells();
    void saveSurfaceCells();
    void saveCellSet(List<label>& listToSave, word fileName);

    // cutting fields in normal direction
    void cutFInBoundaryCells(volVectorField& f);
    void cutUInBoundaryCells(volVectorField& U);
    void cutPhiInBoundaryCells(surfaceScalarField& phi);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "ibInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
