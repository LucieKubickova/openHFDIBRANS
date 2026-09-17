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
    openHFDIBRANS is licensed under the GNU LESSER GENERAL PUBLIC LICENSE
    (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the
    terms and conditions of version 3 of the GNU General Public License,
    supplemented by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIBRANS. If not, see
    <http://www.gnu.org/licenses/lgpl.html>.

Class
    Foam::CLASSNAME

Description

SourceFiles
    CLASSNAMEI.H
    CLASSNAME.C
    CLASSNAMEIO.C

Contributors
	Federico Municchi (2016),
	Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),

\*---------------------------------------------------------------------------*/

#ifndef _Template_H
#define _Template_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class CLASSNAME;
Istream& operator>>(Istream&, CLASSNAME&);
Ostream& operator<<(Ostream&, const CLASSNAME&);


/*---------------------------------------------------------------------------*\
                         Class CLASSNAME Declaration
\*---------------------------------------------------------------------------*/

class CLASSNAME
:
    public baseClassName
{
    // Private Data

        //- Description of data_
        dataType data_;


    // Private Member Functions


public:

    // Static Data Members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        CLASSNAME();

        //- Construct from components
        CLASSNAME(const dataType& data);

        //- Construct from Istream
        CLASSNAME(Istream&);

        //- Copy constructor
        CLASSNAME(const CLASSNAME&);

        //- Disallow default bitwise copy construction
        CLASSNAME(const CLASSNAME&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<CLASSNAME> New();


    //- Destructor
    ~CLASSNAME();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const CLASSNAME&);

        //- Disallow default bitwise assignment
        void operator=(const CLASSNAME&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, CLASSNAME&);
        friend Ostream& operator<<(Ostream&, const CLASSNAME&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CLASSNAMEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
