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

#ifndef ibDirichletBCs_H
#define ibDirichletBCs_H

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
#include "nearWallDist.H"

#include "boundaryCell.H"

using namespace Foam;

//---------------------------------------------------------------------------//
class ibDirichletBCs
{
    private:

    // Private data types ------------------------------------------

    // mesh related parameters
    const fvMesh& mesh_;

    // lambda related things
    const volScalarField& body_;
    scalar thrSurf_;
    bool useYEff_;
    word uTauType_;
    scalar uTauCoeff_;

    // fields
    List<List<scalar>> nutAtIB_;
    List<List<scalar>> uTauAtIB_;

    // references to shared data
    List<DynamicList<boundaryCell>>& boundaryCells_;
    labelField& isBoundaryCell_;

    // dictionaries
    IOdictionary turbulenceProperties_;
    IOdictionary HFDIBDEMDict_;
    IOdictionary fvSchemes_;
    dictionary HFDIBBCsDict_;

    // simulation type
    word simulationType_;

    // saves
    volScalarField yPlusi_;
    volScalarField uTaui_;

    // RAS wall function names
    word nutWF_;
    word kWF_;
    word omegaWF_;
    word epsilonWF_;

    // RAS parameters
    const scalar kappa_;
    const scalar E_;
    const scalar Cmu_;
    const scalar Ceps2_;
    const scalar beta1_;
    scalar Cmu75_;
    scalar Cmu25_;
    scalar Cmu5_;
    scalar yPlusLam_;

    // private functions
    void calcYPlusLam();
    scalar ibZero(scalar dummy)
    {
        return 0.0;
    }

    vector ibZero(vector dummy)
    {
        return vector::zero;
    }

    scalar ibOne(scalar dummy)
    {
        return 1.0;
    }

    vector ibOne(vector dummy)
    {
        return vector::one;
    }

    public:

    // Constructors ----------------------------------------------------
    ibDirichletBCs
    (
        const fvMesh& mesh,
        const volScalarField& body,
        List<DynamicList<boundaryCell>>& boundaryCells,
        labelField& isBoundaryCell
    );

    // Destructors -----------------------------------------------------
    ~ibDirichletBCs();

    // Public methods --------------------------------------------------

    // set size to list
    void setSizeToLists();

    // velocity boundary condition
    void UAtIB(List<vector>& UIB, word BCType);

    // scalar boundary condition
    void TAtIB(List<scalar>& TIB, scalar TIn);

    // calculate friction velocity
    void updateUTauAtIB(volScalarField& k);

    // correct nut at the immersed boundary
    void correctNutAtIB(volScalarField& k, volScalarField& nu);

    // turbulence fields boundary conditions
    void kAtIB
    (
        List<scalar>& kIB,
        volScalarField& k,
        volScalarField& nu
    );
    
    void omegaGAtIB
    (
        List<scalar>& omegaIB,
        List<scalar>& GIB,
        volScalarField::Internal& G,
        const volVectorField& U,
        volScalarField& k,
        volScalarField& nu
    );

    void epsilonGAtIB
    (
        List<scalar>& epsilonIB,
        List<scalar>& GIB,
        volScalarField::Internal& G,
        const volVectorField& U,
        volScalarField& k,
        volScalarField& nu
    );

    // return the immersed yPlus field
    volScalarField& getYPlusi()
    {
        return yPlusi_;
    }

    // return the yPlusLam
    scalar getYPlusLam()
    {
        return yPlusLam_;
    }

    // return E coefficient
    scalar getE()
    {
        return E_;
    }

    // get nut at immersed boundary
    List<List<scalar>>& getNutAtIB()
    {
        return nutAtIB_;
    }

    // post processing of uTau
    void postProcessUTau();

    // save
    void saveUTau();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
