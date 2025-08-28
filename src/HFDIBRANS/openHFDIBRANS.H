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
    scalar thrSurf_;

    // dictionaries
    IOdictionary HFDIBDEMDict_;
    IOdictionary fvSchemes_;
    dictionary HFDIBOuterSchemes_;

    // lists of cells relevant to imposed fields
    List<DynamicList<Tuple3<label,label,label>>> boundaryCells_;
    List<List<Tuple3<scalar,scalar,scalar>>> boundaryDists_;
    List<DynamicList<label>> surfaceCells_;
    List<List<scalar>> surfaceDists_;
    List<DynamicList<label>> internalCells_;
    labelField isBoundaryCell_;

    // interpolation class
    autoPtr<ibInterpolation> ibInterpolation_;

    // wall functions class
    autoPtr<ibDirichletBCs> ibDirichletBCs_;

    // options
    bool save_;
    bool cpDisToInner_;
    bool scaleDisG_;
    scalar scaleCoeff_;
    bool useYEff_;

    // save
    word outDir_;

    public:

    // Constructors ----------------------------------------------------
    openHFDIBRANS(const fvMesh& mesh, const volScalarField& body);

    // Destructors -----------------------------------------------------
    ~openHFDIBRANS();

    // Public methods --------------------------------------------------

    // compute imposed fields
    void computeUi(volVectorField& U, volVectorField& Ui, word surfType);
    void computeKi(volScalarField& k, volScalarField& ki, volScalarField& nu);
    void computeTi(volScalarField& T, volScalarField& Ti, word surfType, scalar TIn);

    // update friction velocity
    void updateUTau(volScalarField& k);
    
    // correct values
    void correctNut(volScalarField& k, volScalarField& nu);
    void correctOmegaG(volScalarField& omega, volScalarField::Internal& G, const volVectorField& U, volScalarField& k, volScalarField& nu, volScalarField& surface);
    void correctEpsilonG(volScalarField& epsilon, volScalarField::Internal& G, const volVectorField& U, volScalarField& k, volScalarField& nu, volScalarField& surface);

    // create and return surface field
    void createBaseSurface(volScalarField& surface, word surfType, scalar boundaryValue);
    void updateSurface(volScalarField& surface, word surfType);

    // correct y for blending functions in kOmegaSST
    void correctY(volScalarField& y);

    // cutting fields in normal direction 
    void cutFInBoundaryCells(volVectorField& f);
    void cutUInBoundaryCells(volVectorField& U);
    void cutPhiInBoundaryCells(surfaceScalarField& phi);

    // enforcing Ui
    void enforceUiInBody(volVectorField& U, volVectorField& Ui);

    // bound fields
    void bound(volScalarField& phi, dimensionedScalar& phiMin);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
