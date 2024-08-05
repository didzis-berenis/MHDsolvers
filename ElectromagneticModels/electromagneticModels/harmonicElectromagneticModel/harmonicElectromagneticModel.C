/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "harmonicElectromagneticModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(harmonicElectromagneticModel, 0);
    defineRunTimeSelectionTable(harmonicElectromagneticModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::harmonicElectromagneticModel::harmonicElectromagneticModel
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    PotEre_(lookupOrConstructScalar(mesh, "PotEre")),
    PotEim_(lookupOrConstructScalar(mesh, "PotEim")),
    
    Jre_(lookupOrConstructVector(mesh, "Jre")),
    Jim_(lookupOrConstructVector(mesh, "Jim")),
    
    Bre_(lookupOrConstructVector(mesh, "Bre")),
    Bim_(lookupOrConstructVector(mesh, "Bim"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::harmonicElectromagneticModel> Foam::harmonicElectromagneticModel::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return electromagneticModel::New<harmonicElectromagneticModel>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::harmonicElectromagneticModel::~harmonicElectromagneticModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::harmonicElectromagneticModel::PotEre()
{
    return PotEre_;
}

Foam::volScalarField& Foam::harmonicElectromagneticModel::PotEim()
{
    return PotEim_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::Jre()
{
    return Jre_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::Jim()
{
    return Jim_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::Bre()
{
    return Bre_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::Bim()
{
    return Bim_;
}


// ************************************************************************* //
