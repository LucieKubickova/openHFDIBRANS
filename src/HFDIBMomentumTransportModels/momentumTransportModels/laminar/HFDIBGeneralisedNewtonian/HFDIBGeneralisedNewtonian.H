/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::HFDIBLaminarModels::HFDIBGeneralisedNewtonian

Description
    Turbulence model for shear-dependent Non-Newtonian flow.

SourceFiles
    HFDIBGeneralisedNewtonian.C

\*---------------------------------------------------------------------------*/

#ifndef HFDIBGeneralisedNewtonian_H
#define HFDIBGeneralisedNewtonian_H

#include "HFDIBLaminarModel.H"
#include "linearViscousStress.H"
#include "generalisedNewtonianViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBLaminarModels
{

/*---------------------------------------------------------------------------*\
                    Class HFDIBGeneralisedNewtonian Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class HFDIBGeneralisedNewtonian
:
    public linearViscousStress<HFDIBLaminarModel<BasicMomentumTransportModel>>
{

protected:

    // Protected data

        //- Run-time selectable generalised Newtonian viscosity model
        autoPtr<Foam::laminarModels::generalisedNewtonianViscosityModel> viscosityModel_;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("generalizedNewtonian"); // Note: do not use HFDIBGeneralisedNewtonian


    // Constructors

        //- Construct from components
        HFDIBGeneralisedNewtonian
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<HFDIBGeneralisedNewtonian> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~HFDIBGeneralisedNewtonian()
    {}


    // Member Functions

        //- Read momentumTransport dictionary
        virtual bool read();

        //- Return the effective viscosity
        //  i.e. the HFDIBGeneralisedNewtonian viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const;

        //- Correct the HFDIBGeneralisedNewtonian viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBLaminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBGeneralisedNewtonian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
