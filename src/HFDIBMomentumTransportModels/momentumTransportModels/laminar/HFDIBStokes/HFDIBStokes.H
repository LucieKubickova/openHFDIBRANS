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

Class
    Foam::HFDIBLaminarModels::HFDIBStokes

Description
    Turbulence model for HFDIBStokes flow.

SourceFiles
    HFDIBStokes.C

\*---------------------------------------------------------------------------*/

#ifndef HFDIBStokes_H
#define HFDIBStokes_H

#include "HFDIBLaminarModel.H"
#include "linearViscousStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBLaminarModels
{

/*---------------------------------------------------------------------------* \
                           Class HFDIBStokes Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class HFDIBStokes
:
    public linearViscousStress<HFDIBLaminarModel<BasicMomentumTransportModel>>
{

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("Stokes"); // Note: do not use HFDIBStokes


    // Constructors

        //- Construct from components
        HFDIBStokes
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<HFDIBStokes> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );


    //- Destructor
    virtual ~HFDIBStokes()
    {}


    // Member Functions

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const;

        //- Read momentumTransport dictionary
        virtual bool read();

        //- Return the turbulence viscosity, i.e. 0 for HFDIBStokes flow
        virtual tmp<volScalarField> nut() const;

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const;

        //- Return the effective viscosity, i.e. the HFDIBStokes viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const;

        //- Return the turbulence kinetic energy, i.e. 0 for HFDIBStokes flow
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate,
        //  i.e. 0 for HFDIBStokes flow
        virtual tmp<volScalarField> epsilon() const;

        //- Return the stress tensor [m^2/s^2], i.e. 0 for HFDIBStokes flow
        virtual tmp<volSymmTensorField> sigma() const;

        //- Correct the HFDIBStokes viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBLaminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBStokes.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
