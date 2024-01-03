/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::HFDIBLaminarModel

Description
    Templated abstract base class for HFDIBLaminar transport models

SourceFiles
    HFDIBLaminarModel.C

\*---------------------------------------------------------------------------*/

#ifndef HFDIBLaminarModel_H
#define HFDIBLaminarModel_H

#include "HFDIBmomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HFDIBLaminarModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicHFDIBMomentumTransportModel>
class HFDIBLaminarModel
:
    public BasicHFDIBMomentumTransportModel
{

protected:

    // Protected data

        //- HFDIBLaminar coefficients dictionary
        dictionary HFDIBLaminarDict_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs(const word& type);


public:

    typedef typename BasicHFDIBMomentumTransportModel::alphaField alphaField;
    typedef typename BasicHFDIBMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("laminar"); // Note: do not change to HFDIBLaminar


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            HFDIBLaminarModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const viscosity& viscosity
            ),
            (alpha, rho, U, alphaRhoPhi, phi, viscosity)
        );


    // Constructors

        //- Construct from components
        HFDIBLaminarModel
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );

        //- Disallow default bitwise copy construction
        HFDIBLaminarModel(const HFDIBLaminarModel&) = delete;


    // Selectors

        //- Return a reference to the selected HFDIBLaminar model
        static autoPtr<HFDIBLaminarModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~HFDIBLaminarModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Return the turbulence viscosity, i.e. 0 for HFDIBLaminar flow
        virtual tmp<volScalarField> nut() const;

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const;

        //- Return the effective viscosity, i.e. the HFDIBLaminar viscosity
        virtual tmp<volScalarField> nuEff() const = 0;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const = 0;

        //- Return the turbulence kinetic energy, i.e. 0 for HFDIBLaminar flow
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate,
        //  i.e. 0 for HFDIBLaminar flow
        virtual tmp<volScalarField> epsilon() const;

        //- Return the turbulence specific dissipation rate,
        //  i.e. 0 for laminar flow
        virtual tmp<volScalarField> omega() const;

        //- Return the stress tensor [m^2/s^2], i.e. 0 for HFDIBLaminar flow
        virtual tmp<volSymmTensorField> sigma() const;

        //- Correct the HFDIBLaminar viscosity
        virtual void correct();

        //- second image of the correct function to accept arguments
        virtual void correct(openHFDIBRANS& HFDIBRANS);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HFDIBLaminarModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBLaminarModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
