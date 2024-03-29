/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::HFDIBLaminarModels::HFDIBPTT

Description
    HFDIBPTT model for viscoelasticity using the upper-convected time
    derivative of the stress tensor with support for multiple modes.

    Reference:
    \verbatim
        Thien, N. P., & Tanner, R. I. (1977).
        A new constitutive equation derived from network theory.
        Journal of Non-Newtonian Fluid Mechanics, 2(4), 353-365.
    \endverbatim

See also
    Foam::HFDIBLaminarModels::HFDIBMaxwell
    Foam::HFDIBLaminarModels::Giesekus

SourceFiles
    HFDIBPTT.C

\*---------------------------------------------------------------------------*/

#ifndef HFDIBPTT_H
#define HFDIBPTT_H

#include "HFDIBMaxwell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HFDIBLaminarModels
{

/*---------------------------------------------------------------------------*\
                           Class HFDIBPTT Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class HFDIBPTT
:
    public HFDIBMaxwell<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Mode coefficients

            PtrList<dimensionedScalar> epsilons_;


    // Protected Member Functions

        virtual tmp<fvSymmTensorMatrix> sigmaSource
        (
            const label modei,
            volSymmTensorField& sigma
        ) const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("PTT"); // Note: do not use HFDIBPTT


    // Constructors

        //- Construct from components
        HFDIBPTT
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        HFDIBPTT(const HFDIBPTT&) = delete;



    //- Destructor
    virtual ~HFDIBPTT()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HFDIBPTT&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HFDIBLaminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HFDIBPTT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
