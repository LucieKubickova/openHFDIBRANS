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

    openHFDIBRANS-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

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
    extended by connection with (level-set) DEM (Isoz and Sourek, 2020)

SourceFiles
    openHFDIBRANSDEM.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef openHFDIBRANS_H
#define openHFDIBRANS_H

#include "ibInterpolation.H"
#include "ibDirichletBCs.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
class openHFDIBRANS
{
    private:

    // Private data types ------------------------------------------

    // mesh related parameters
    const fvMesh& mesh_;

    // lambda related
    const volScalarField& body_;
    volScalarField surface_;

    // turbulence related
    word simulationType_;

    // dictionaries
    IOdictionary HFDIBDEMDict_;

    // lists of cells relevant to imposed fields
    DynamicList<Tuple2<label,label>> boundaryCells_;
    List<Tuple2<scalar,scalar>> boundaryDists_;

    // interpolation class
    ibInterpolation ibInterpolation_;

    // wall functions class
    ibDirichletBCs ibDirichletBCs_;

    // save
    bool save_;
    word outDir_;

    public:

    // Constructors ----------------------------------------------------
    openHFDIBRANS(const fvMesh& mesh, const volScalarField& body, word simulationType);

    // Destructors -----------------------------------------------------
    ~openHFDIBRANS();

    // Public methods --------------------------------------------------
    
    // compute imposed fields
    void computeUi(volVectorField& U, volVectorField& Ui, const volScalarField& k, const volScalarField& nu);
    void computeKi(volScalarField& k, volScalarField& ki, volScalarField& nu);

    // correct fields
    void correctOmegaG(volScalarField& omega, volScalarField::Internal& G, const volVectorField& U, volScalarField& k, volScalarField& nu);

    // return the surface field
    volScalarField& surface()
    {
        return surface_;
    }

    // return the surface threshold
    scalar getSurfaceThreshold()
    {
        return ibInterpolation_.getSurfaceThreshold();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
