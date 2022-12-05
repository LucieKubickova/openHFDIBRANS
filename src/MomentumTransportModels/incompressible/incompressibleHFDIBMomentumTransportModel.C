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

#include "incompressibleHFDIBMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleHFDIBMomentumTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleHFDIBMomentumTransportModel::incompressibleHFDIBMomentumTransportModel
(
    const geometricOneField&,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi
)
:
    HFDIBmomentumTransportModel
    (
        U,
        alphaRhoPhi,
        phi
    )
{}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleHFDIBMomentumTransportModel::mu() const
{
    return nu();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleHFDIBMomentumTransportModel::mu(const label patchi) const
{
    return nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleHFDIBMomentumTransportModel::mut() const
{
    return nut();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleHFDIBMomentumTransportModel::mut(const label patchi) const
{
    return nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleHFDIBMomentumTransportModel::muEff() const
{
    return nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleHFDIBMomentumTransportModel::muEff(const label patchi) const
{
    return nuEff(patchi);
}


// ************************************************************************* //
