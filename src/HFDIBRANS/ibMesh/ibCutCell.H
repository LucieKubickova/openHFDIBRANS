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
    cutCell class to be used for hfdibRANS turbulence models

Contributors
    Martin Isoz (2019-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef ibCutCell_H
#define ibCutCell_H

#include "fvMesh.H"
#include "List.H"

namespace Foam
{
//---------------------------------------------------------------------------//
class ibCutCell
{
    private:

    // Private data types ------------------------------------------
    const fvMesh& mesh_;
    
    const vector& surfNorm_;
    const point& surfPoint_;
    const cell& inCell_;
    
    
    //~ DynamicList<point> v_;
    DynamicList<DynamicList<point>> v_;
    
    
    //~ faceList faces_;
    //~ DynamicList<List<label>> faces_;                                           //face class is too complicated
    DynamicList<DynamicList<point>> faces_;                                           //face class is too complicated
    List<vector> faceCentres_;
    List<vector> faceAreas_;
    
    vector C_;                                                          //cutCell centre
    scalar V_;                                                          //cutCell volume
    vector CIn_;                                                        //input cell centre
    scalar VIn_;                                                        //input cell volume
        
    label type_;                                                        //-1 .. in, 1 .. out, 0 .. cut
    label nCutFaces_;

    // private member functions
    void makeFaces(const List<point>& vertices);
    
    void makeCellCentreAndVol();
    void makeFaceCentresAndAreas();
    
    void syncNormals();
    
    protected:

    public:

    // Constructors ----------------------------------------------------
    ibCutCell(
        const fvMesh& mesh,                                             //reference to fvMesh
        const vector& surfNorm,                                         //normal to the cutting plane
        const point& surfPoint,                                         //point on the cutting plane
        const cell& inCell                                              //input cell
    );

    // Destructors -----------------------------------------------------
    ~ibCutCell();

    // Public methods --------------------------------------------------
    
    // getters
    vector C(){return C_;};                                             //return cutCell centre
    scalar V(){return V_;};                                             //return cutCell volume
    vector S(){return surfPoint_;};                                     //return surfPoint
    
    vector CIn(){return CIn_;};                                         //return input cell centre
    scalar VIn(){return VIn_;};                                         //return input cell volume
    
    List<vector> Sf(){return faceAreas_;};
    List<vector> n();
    vector Ss(){return surfNorm_;};                                     //return cutting surface normal
    
    DynamicList<DynamicList<point>> faces(){return faces_;};            //return cutCell faces
    
    label status(){return type_;};
    
    scalar yOrtho();                                                    //returns orthogonal distance from C() to the wall
    
    scalar yOrthoEst();                                                 //returns 1/2 of distance from S() to the oposite face
    
    //~ faceList faces();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
