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
    DynamicList<Tuple2<label,label>> boundaryCells_;
    List<Tuple2<scalar,scalar>> boundaryDists_;
    labelField isBoundaryCell_;

    // interpolation class
    ibInterpolation ibInterpolation_;

    // wall functions class
    ibDirichletBCs ibDirichletBCs_;

    // options
    bool save_;
    bool cpOmegaToInner_;
    bool scaleDisG_;

    // save
    word outDir_;

    public:

    // Constructors ----------------------------------------------------
    openHFDIBRANS(const fvMesh& mesh, const volScalarField& body);

    // Destructors -----------------------------------------------------
    ~openHFDIBRANS();

    // Public methods --------------------------------------------------

    // compute imposed fields
    void computeUi(volVectorField& U, volVectorField& Ui);
    void computeKi(volScalarField& k, volScalarField& ki, volScalarField& nu);

    // update friction velocity
    void updateUTau(volScalarField& k);
    
    // correct values
    void correctNut(volScalarField& k, volScalarField& nu);
    void correctOmegaG(volScalarField& omega, volScalarField::Internal& G, const volVectorField& U, volScalarField& k, volScalarField& nu, volScalarField& surface);
    void correctEpsilonG(volScalarField& epsilon, volScalarField::Internal& G, const volVectorField& U, volScalarField& k, volScalarField& nu, volScalarField& surface);

    // create and return surface field
    void createBaseSurface(volScalarField& surface, word surfType, scalar boundaryValue);
    void updateSurface(volScalarField& surface, word surfType);

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
