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

using namespace Foam;

//---------------------------------------------------------------------------//
class ibDirichletBCs
{
    private:

    // Private data types ------------------------------------------

    // mesh related parameters
    const fvMesh& mesh_;

    // references to shared data
    DynamicList<Tuple2<label,label>>& boundaryCells_;
    List<Tuple2<scalar,scalar>>& boundaryDists_;
    DynamicList<label>& boundaryFaces_;
    List<bool>& isWallCell_;

    // dictionaries
    IOdictionary HFDIBDEMDict_;
    dictionary HFDIBBCsDict_;

    // simulation type
    word simulationType_;

    // saves
    List<scalar> yPlusi_;
    List<scalar> ULogScales_;
    List<scalar> kLogScales_;

    // RAS wall function names
    word kWF_;
    word omegaWF_;

    // RAS parameters
    const scalar kappa_;
    const scalar E_;
    const scalar Cmu_;
    const scalar Ceps2_;
    const scalar beta1_;
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
    ibDirichletBCs(const fvMesh& mesh, word simulationType, DynamicList<Tuple2<label,label>>& boundaryCells, List<Tuple2<scalar,scalar>>& boundaryDists, DynamicList<label>& boundaryFaces, List<bool>& isWallCell);

    // Destructors -----------------------------------------------------
    ~ibDirichletBCs();

    // Public methods --------------------------------------------------
    
    // velocity boudnary condition
    template <typename Type>
    void UAtIB(List<Type>& UIB, word BCType, const volScalarField& k, const volScalarField& nu);

    // pressure boundary condition
    void pAtIB(List<scalar>& pIB, word BCType, const volScalarField& p);

    // turbulence fields boundary conditions
    void omegaGAtIB(List<scalar>& omegaIB, List<scalar>& GIB, volScalarField::Internal& G, const volVectorField& U, volScalarField& k, volScalarField& nu);
    void kAtIB(List<scalar>& kIB, volScalarField& k, volScalarField& nu);

    // flux boundary condition
    void phiAtIB(List<scalar>& phiIB, word BCType);

    // return the log scales for velocity
    List<scalar>& getULogScales()
    {
        return ULogScales_;
    }

    // return the log scales for turbulent kinetic energy
    List<scalar>& getkLogScales()
    {
        return kLogScales_;
    }

    // return the immersed yPlus field
    List<scalar>& getYPlusi()
    {
        return yPlusi_;
    }

    // return the yPlusLam
    scalar getYPlusLam()
    {
        return yPlusLam_;
    }

    // correct lists
    void correctLists()
    {
        yPlusi_.setSize(boundaryCells_.size());
        ULogScales_.setSize(boundaryCells_.size());
        kLogScales_.setSize(boundaryCells_.size());
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "ibDirichletBCsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
