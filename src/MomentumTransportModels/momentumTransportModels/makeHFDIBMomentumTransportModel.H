/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#define makeHFDIBMomentumTransportModelTypes(                                  \
    Alpha, Rho, baseModel, BaseModel, Transport)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel<Transport> Transport##BaseModel;                     \
        typedef HFDIBRASModel<Transport##BaseModel> HFDIBRAS##Transport##BaseModel;\
    }


#define makeBaseHFDIBMomentumTransportModel(                                   \
    Alpha, Rho, baseModel, BaseModel, Transport)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef HFDIBMomentumTransportModel                                         \
        <                                                                      \
            Alpha,                                                             \
            Rho,                                                               \
            baseModel,                                                         \
            Transport                                                          \
        > Transport##baseModel;                                                \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            Transport##baseModel,                                              \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(HFDIBRAS##Transport##BaseModel, 0);\
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (HFDIBRAS##Transport##BaseModel, dictionary);                          \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##baseModel,                                              \
            HFDIBRAS##Transport##BaseModel,                                    \
            dictionary                                                         \
        );                                                                     \
    }

#define makeTemplatedHFDIBMomentumTransportModel(BaseModel, SType, Type)       \
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


#define makeHFDIBMomentumTransportModel(BaseModel, SType, Type)                \
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
