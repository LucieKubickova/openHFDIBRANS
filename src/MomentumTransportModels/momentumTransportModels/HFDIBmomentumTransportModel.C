/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "HFDIBmomentumTransportModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HFDIBmomentumTransportModel, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HFDIBmomentumTransportModel::HFDIBmomentumTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi
)
: momentumTransportModel(U, alphaRhoPhi, phi)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::HFDIBmomentumTransportModel::phi() const
{
    return phi_;
}


bool Foam::HFDIBmomentumTransportModel::read()
{
    return regIOobject::read();
}


void Foam::HFDIBmomentumTransportModel::validate()
{}


void Foam::HFDIBmomentumTransportModel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


void Foam::HFDIBmomentumTransportModel::correct(openHFDIBRANS& HFDIBRANS)
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


// ************************************************************************* //
