/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMesh.H"

//~ template<class HFDIBMomentumTransportModel>
//~ inline Foam::autoPtr<HFDIBMomentumTransportModel> Foam::HFDIBmomentumTransportModel::New
//~ (
    //~ const typename HFDIBMomentumTransportModel::alphaField& alpha,
    //~ const typename HFDIBMomentumTransportModel::rhoField& rho,
    //~ const volVectorField& U,
    //~ const surfaceScalarField& alphaRhoPhi,
    //~ const surfaceScalarField& phi,
    //~ const viscosity& viscosity
//~ )
//~ {
    //~ const word modelType
    //~ (
        //~ IOdictionary
        //~ (
            //~ HFDIBmomentumTransportModel::readModelDict
            //~ (
                //~ U.db(),
                //~ alphaRhoPhi.group()
            //~ )
        //~ ).lookup("simulationType")
    //~ );

    //~ Info<< "Selecting turbulence model type " << modelType << endl;

    //~ typename HFDIBMomentumTransportModel::dictionaryConstructorTable::iterator
        //~ cstrIter =
        //~ HFDIBMomentumTransportModel::dictionaryConstructorTablePtr_->find(modelType);

    //~ if
    //~ (
        //~ cstrIter
     //~ == HFDIBMomentumTransportModel::dictionaryConstructorTablePtr_->end()
    //~ )
    //~ {
        //~ FatalErrorInFunction
            //~ << "Unknown " << HFDIBMomentumTransportModel::typeName << " type "
            //~ << modelType << nl << nl
            //~ << "Valid " << HFDIBMomentumTransportModel::typeName << " types:" << endl
            //~ << HFDIBMomentumTransportModel::dictionaryConstructorTablePtr_
               //~ ->sortedToc()
            //~ << exit(FatalError);
    //~ }

    //~ return autoPtr<HFDIBMomentumTransportModel>
    //~ (
        //~ cstrIter()(alpha, rho, U, alphaRhoPhi, phi, viscosity)
    //~ );
//~ }


// ************************************************************************* //
