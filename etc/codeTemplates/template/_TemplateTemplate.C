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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*), Lucie Kubíčková (2021-*),

\*---------------------------------------------------------------------------*/

#include "CLASSNAME.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<TemplateClassArgument>
const dataType CLASSNAME<TemplateArgument>::staticData();


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<TemplateClassArgument>
CLASSNAME<TemplateArgument>::CLASSNAME()
:
    baseClassName(),
    data_()
{}


template<TemplateClassArgument>
CLASSNAME<TemplateArgument>::CLASSNAME(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


template<TemplateClassArgument>
CLASSNAME<TemplateArgument>::CLASSNAME
(
    const CLASSNAME<TemplateArgument>&
)
:
    baseCLASSNAME(),
    data_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<TemplateClassArgument>
autoPtr<CLASSNAME<TemplateArgument>>
CLASSNAME<TemplateArgument>::New()
{
    return autoPtr<CLASSNAME<TemplateArgument>>
    (
        new CLASSNAME<TemplateArgument>
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<TemplateClassArgument>
CLASSNAME<TemplateArgument>::~CLASSNAME()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<TemplateClassArgument>
void CLASSNAME<TemplateArgument>::operator=
(
    const CLASSNAME<TemplateArgument>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
