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

#include "IncompressibleHFDIBMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>::
IncompressibleHFDIBMomentumTransportModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport
)
:
    HFDIBMomentumTransportModel
    <
        geometricOneField,
        geometricOneField,
        incompressibleHFDIBMomentumTransportModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>>
Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const TransportModel& transport
)
{
    return autoPtr<IncompressibleHFDIBMomentumTransportModel>
    (
        static_cast<IncompressibleHFDIBMomentumTransportModel*>(
        HFDIBMomentumTransportModel
        <
            geometricOneField,
            geometricOneField,
            incompressibleHFDIBMomentumTransportModel,
            TransportModel
        >::New
        (
            geometricOneField(),
            geometricOneField(),
            U,
            phi,
            phi,
            transport
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>::devSigma() const
{
    return devTau();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>::divDevSigma
(
    volVectorField& U
) const
{
    return divDevTau(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>::
devTau() const
{
    NotImplemented;

    return devSigma();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>::
divDevTau
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleHFDIBMomentumTransportModel<TransportModel>::
divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


// ************************************************************************* //
