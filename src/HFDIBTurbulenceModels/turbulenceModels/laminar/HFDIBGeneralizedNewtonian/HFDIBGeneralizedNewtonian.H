/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Class
    Foam::HFDIBLaminarModels::HFDIBGeneralizedNewtonian

Description
    Turbulence model for shear-dependent Non-Newtonian flow.

SourceFiles
    HFDIBGeneralizedNewtonian.C

\*---------------------------------------------------------------------------*/

#ifndef HFDIBGeneralizedNewtonian_H
#define HFDIBGeneralizedNewtonian_H

#include "HFDIBLaminarModel.H"
#include "linearViscousStress.H"
#include "generalizedNewtonianViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBLaminarModels
{

/*---------------------------------------------------------------------------*\
                    Class HFDIBGeneralizedNewtonian Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class HFDIBGeneralizedNewtonian
:
    public linearViscousStress<HFDIBLaminarModel<BasicMomentumTransportModel>>
{
protected:

    // Protected Data

        //- Run-time selectable non-Newtonian viscosity model
        autoPtr<Foam::laminarModels::generalizedNewtonianViscosityModel> viscosityModel_;

        //- The non-Newtonian viscosity field
        volScalarField nu_;


    // Protected Member Functions

        virtual tmp<volScalarField> strainRate() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("HFDIBGeneralizedNewtonian");


    // Constructors

        //- Construct from components
        HFDIBGeneralizedNewtonian
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<HFDIBGeneralizedNewtonian> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    //- Destructor
    virtual ~HFDIBGeneralizedNewtonian() = default;


    // Member Functions

        //- Read turbulence (momentumTransport) dictionary
        virtual bool read();

        //- Return the turbulence viscosity,
        virtual tmp<volScalarField> nut() const;

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const;

        //- Correct the HFDIBGeneralizedNewtonian viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBLaminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBGeneralizedNewtonian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
