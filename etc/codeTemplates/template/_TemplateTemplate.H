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

#ifndef _TemplateTemplate_H
#define _TemplateTemplate_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class someClass;


// Forward declaration of friend functions and operators

template<TemplateClassArgument>
class CLASSNAME;

template<TemplateClassArgument>
Istream& operator>>(Istream&, CLASSNAME<TemplateArgument>&);

template<TemplateClassArgument>
Ostream& operator<<(Ostream&, const CLASSNAME<TemplateArgument>&);


/*---------------------------------------------------------------------------*\
                         Class CLASSNAME Declaration
\*---------------------------------------------------------------------------*/

template<TemplateClassArgument>
class CLASSNAME
:
    public baseClassName
{
    // Private Data

        dataType data_;


    // Private Member Functions


public:

    // Static Data Members

        //- Static data someStaticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        CLASSNAME();

        //- Construct from components
        CLASSNAME(const dataType& data);

        //- Construct from Istream
        CLASSNAME(Istream&);

        //- Copy constructor
        CLASSNAME(const CLASSNAME<TemplateArgument>&);

        //- Disallow default bitwise copy construction
        CLASSNAME(const CLASSNAME<TemplateArgument>&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<CLASSNAME<TemplateArgument>> New();


    //- Destructor
    ~CLASSNAME();


    // Member Functions

    // Member Operators

        void operator=(const CLASSNAME<TemplateArgument>&);

        //- Disallow default bitwise assignment
        void operator=(const CLASSNAME<TemplateArgument>&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>> <TemplateArgument>
        (Istream&, CLASSNAME<TemplateArgument>&);

        friend Ostream& operator<< <TemplateArgument>
        (Ostream&, const CLASSNAME<TemplateArgument>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CLASSNAMEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CLASSNAME.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
