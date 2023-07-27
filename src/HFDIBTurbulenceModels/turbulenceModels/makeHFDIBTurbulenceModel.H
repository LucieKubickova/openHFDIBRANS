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

//- Define turbulence model typedefs, with EddyDiffusivity
#define defineHFDIBTurbulenceModelTypes(                                            \
    Alpha, Rho, baseModel, BaseModel, Transport)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef HFDIBTurbulenceModel                                                \
        <                                                                      \
            Alpha,                                                             \
            Rho,                                                               \
            baseModel,                                                         \
            Transport                                                          \
        > Transport##baseModel;                                                \
                                                                               \
        typedef BaseModel<Transport>                                           \
            Transport##BaseModel;                                              \
                                                                               \
        typedef HFDIBLaminarModel<Transport##BaseModel>                             \
            HFDIBLaminar##Transport##BaseModel;                                     \
                                                                               \
        typedef HFDIBRASModel<Transport##BaseModel>                                 \
            HFDIBRAS##Transport##BaseModel;                                         \
    }

// Legacy name
#define makeHFDIBTurbulenceModelTypes(a, b, c, d, e) \
    defineHFDIBTurbulenceModelTypes(a, b, c, d, e)


#define makeBaseHFDIBTurbulenceModel(Alpha, Rho, baseModel, BaseModel, Transport)   \
    /* HFDIBTurbulence typedefs */                                                  \
    defineHFDIBTurbulenceModelTypes(Alpha, Rho, baseModel, BaseModel, Transport)    \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        /* HFDIBTurbulence selection table */                                       \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            Transport##baseModel,                                              \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        /* Laminar models */                                                   \
        defineNamedTemplateTypeNameAndDebug(HFDIBLaminar##Transport##BaseModel, 0); \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (HFDIBLaminar##Transport##BaseModel, dictionary);                           \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##baseModel,                                              \
            HFDIBLaminar##Transport##BaseModel,                                     \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        /* HFDIBRAS models */                                                       \
        defineNamedTemplateTypeNameAndDebug(HFDIBRAS##Transport##BaseModel, 0);     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (HFDIBRAS##Transport##BaseModel, dictionary);                               \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##baseModel,                                              \
            HFDIBRAS##Transport##BaseModel,                                         \
            dictionary                                                         \
        );                                                                     \
    }


#define makeTemplatedHFDIBTurbulenceModel(BaseModel, SType, Type)                   \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SType##Models::Type<Foam::BaseModel>, 0);                       \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##Models                                                \
        {                                                                      \
            typedef Type<BaseModel> Type##SType##BaseModel;                    \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##BaseModel,                                              \
                Type##SType##BaseModel,                                        \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


#define makeHFDIBTurbulenceModel(BaseModel, SType, Type)                            \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##Models                                                \
        {                                                                      \
            defineTypeNameAndDebug(Type, 0);                                   \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##BaseModel,                                              \
                Type,                                                          \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


// ************************************************************************* //
