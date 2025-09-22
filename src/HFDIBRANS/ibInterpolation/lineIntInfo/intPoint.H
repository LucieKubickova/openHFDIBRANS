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
    struct for interpolation point

SourceFiles

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2021-*)
\*---------------------------------------------------------------------------*/

#ifndef intPoint_H
#define intPoint_H

namespace Foam
{
    struct intPoint
    {
        // interpolation point
        point iPoint_ = vector::zero;

        // cell label
        label iCell_ = 0;

        // current processor
        label iProc_ = -1;

        // processor of origin
        label oProc_ = -1;

        // Constructor
        intPoint():
        iPoint_(vector::zero),
        iCell_(0),
        iProc_(-1),
        oProc_(-1)
        {};

        // Copy constructor
        intPoint
        (
            point iPoint,
            label iCell,
            label iProc,
            label oProc
        ):
        iPoint_(iPoint),
        iCell_(iCell),
        iProc_(iProc),
        oProc_(oProc)
        {};

        // Destructor
        ~intPoint(){};

        // Operators
        friend bool operator==(const intPoint& aP, const intPoint& bP)
        {
            return aP.iPoint_ == bP.iPoint_ &&
                   aP.iCell_ == bP.iCell_ &&
                   aP.iProc_ == bP.iProc_ &&
                   aP.oProc_ == bP.oProc_;
        }

        friend bool operator!=(const intPoint& aP, const intPoint& bP)
        {
            return !(aP==bP);
        }

        friend Foam::Ostream& operator<<(Foam::Ostream& os, const intPoint& iP)
        {
            return os << iP.iPoint_ << token::SPACE
                      << iP.iCell_ << token::SPACE
                      << iP.iProc_ << token::SPACE
                      << iP.oProc_;
        }
        
        friend Foam::Istream& operator>>(Foam::Istream& is, intPoint& iP)
        {
            return is >> iP.iPoint_ >> iP.iCell_ >> iP.iProc_ >> iP.oProc_;
        }
    };
}

#endif
