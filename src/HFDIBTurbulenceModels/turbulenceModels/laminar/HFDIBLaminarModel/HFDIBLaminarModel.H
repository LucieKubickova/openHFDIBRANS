/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "HFDIBTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HFDIBLaminarModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicHFDIBTurbulenceModel>
class HFDIBLaminarModel
:
    public BasicHFDIBTurbulenceModel
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

        //- No copy construct
        HFDIBLaminarModel(const HFDIBLaminarModel&) = delete;

        //- No copy assignment
        void operator=(const HFDIBLaminarModel&) = delete;


public:

    typedef typename BasicHFDIBTurbulenceModel::alphaField alphaField;
    typedef typename BasicHFDIBTurbulenceModel::rhoField rhoField;
    typedef typename BasicHFDIBTurbulenceModel::transportModel transportModel;


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
                const transportModel& transport,
                const word& propertiesName
            ),
            (alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
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
            const transportModel& transport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected HFDIBLaminar model
        static autoPtr<HFDIBLaminarModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = HFDIBturbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~HFDIBLaminarModel() = default;


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();


        // Access

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
            virtual tmp<volScalarField> nuEff() const;

            //- Return the effective viscosity on patch
            virtual tmp<scalarField> nuEff(const label patchi) const;

            //- Return the turbulence kinetic energy, i.e. 0 for HFDIBLaminar flow
            virtual tmp<volScalarField> k() const;

            //- Return the turbulence kinetic energy dissipation rate,
            //  i.e. 0 for HFDIBLaminar flow
            virtual tmp<volScalarField> epsilon() const;

            //- Return the specific dissipation rate, i.e. 0 for HFDIBLaminar flow
            virtual tmp<volScalarField> omega() const;

            //- Return the Reynolds stress tensor, i.e. 0 for HFDIBLaminar flow
            virtual tmp<volSymmTensorField> R() const;

            //- Correct the HFDIBLaminar transport
            virtual void correct();

            //- second image of the correct function to accept arguments
            virtual void correct(openHFDIBRANS& HFDIBRANS);
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
