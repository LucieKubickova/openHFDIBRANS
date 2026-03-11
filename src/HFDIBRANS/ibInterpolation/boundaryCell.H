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
    struct for boundary cell

SourceFiles

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef boundaryCell_H
#define boundaryCell_H

namespace Foam
{
    struct boundaryCell
    {
        // cell labels
        label bCell_ = 0; // boundary cell
        label iCell_ = 0; // inner cell
        label fCell1_ = 0; // first free stream cell
        label fCell2_ = 0; // second free stream cell

        // proc labels
        label iProc_ = -1; // proc of the inner cell
        label iFace_ = -1; // label of the shared face at the proc-proc boundary
        label fProc1_ = -1; // proc of the first free stream cell
        label fProc2_ = -1; // proc of the second free stream cell

        // interpolation info
        vector sNorm_ = vector::zero;
        point sPoint_ = vector::zero;
        point fPoint1_ = vector::zero;
        point fPoint2_ = vector::zero;

        // distances
        scalar sigma_ = 0.0; // signed distance 
        scalar yOrtho_ = 0.0; // orthogonal distance
        scalar yEff_ = 0.0; // effective distance

        // Constructor
        boundaryCell():
        bCell_(0),
        iCell_(0),
        fCell1_(0),
        fCell2_(0),
        iProc_(-1),
        iFace_(-1),
        fProc1_(-1),
        fProc2_(-1),
        sNorm_(vector::zero),
        sPoint_(vector::zero),
        fPoint1_(vector::zero),
        fPoint2_(vector::zero),
        sigma_(0.0),
        yOrtho_(0.0),
        yEff_(0.0)
        {};

        // Copy constructor
        boundaryCell
        (
            label bCell,
            label iCell,
            label fCell1,
            label fCell2,
            label iProc,
            label iFace,
            label fProc1,
            label fProc2,
            vector sNorm,
            point sPoint,
            point fPoint1,
            point fPoint2,
            scalar sigma,
            scalar yOrtho,
            scalar yEff
        ):
        bCell_(bCell),
        iCell_(iCell),
        fCell1_(fCell1),
        fCell2_(fCell2),
        iProc_(iProc),
        iFace_(iFace),
        fProc1_(fProc1),
        fProc2_(fProc2),
        sNorm_(sNorm),
        fPoint1_(fPoint1),
        fPoint2_(fPoint2),
        sigma_(sigma),
        yOrtho_(yOrtho),
        yEff_(yEff)
        {};

        // Destructor
        ~boundaryCell(){};
    };
}
#endif
