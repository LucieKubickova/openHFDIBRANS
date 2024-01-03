/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "makeIncompressibleHFDIBMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeBaseHFDIBMomentumTransportModel
(
    geometricOneField,
    geometricOneField,
    incompressibleHFDIBMomentumTransportModel
);


// -------------------------------------------------------------------------- //
// HFDIBLaminar models
// -------------------------------------------------------------------------- //

#include "../../momentumTransportModels/lnInclude/HFDIBStokes.H"
makeHFDIBLaminarModel(HFDIBStokes);

#include "../../momentumTransportModels/lnInclude/HFDIBGeneralisedNewtonian.H"
makeHFDIBLaminarModel(HFDIBGeneralisedNewtonian);

#include "../../momentumTransportModels/lnInclude/HFDIBMaxwell.H"
makeHFDIBLaminarModel(HFDIBMaxwell);

#include "../../momentumTransportModels/lnInclude/HFDIBGiesekus.H"
makeHFDIBLaminarModel(HFDIBGiesekus);

#include "../../momentumTransportModels/lnInclude/HFDIBPTT.H"
makeHFDIBLaminarModel(HFDIBPTT);


// -------------------------------------------------------------------------- //
// HFDIBRAS models
// -------------------------------------------------------------------------- //

#include "../../momentumTransportModels/lnInclude/HFDIBKOmega.H"
makeHFDIBRASModel(HFDIBKOmega);

#include "../../momentumTransportModels/lnInclude/HFDIBKOmegaSST.H"
makeHFDIBRASModel(HFDIBKOmegaSST);

#include "../../momentumTransportModels/lnInclude/HFDIBKEpsilon.H"
makeHFDIBRASModel(HFDIBKEpsilon);

#include "../../momentumTransportModels/lnInclude/HFDIBRealizableKE.H"
makeHFDIBRASModel(HFDIBRealizableKE);

// ************************************************************************* //
