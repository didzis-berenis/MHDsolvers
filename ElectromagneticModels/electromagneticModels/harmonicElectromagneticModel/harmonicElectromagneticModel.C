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
    // Add to runtime selection Table
    /*addToRunTimeSelectionTable
    (
        electromagneticModel,
        harmonicElectromagneticModel,
        fvMesh
    );*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::harmonicElectromagneticModel::harmonicElectromagneticModel
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    electromagneticModel(mesh),
    PotEre_(lookupOrConstructScalar(mesh, "PotEre")),
    PotEim_(lookupOrConstructScalar(mesh, "PotEim")),
    Jre_(lookupOrConstructVector(mesh, "Jre")),
    Jim_(lookupOrConstructVector(mesh, "Jim")),
    Bre_(lookupOrConstructVector(mesh, "Bre")),
    Bim_(lookupOrConstructVector(mesh, "Bim")),
    //Get boundary conditions from J
    deltaJre_
    (
        lookupOrConstructVector
        (
            mesh,
            "deltaJre",
            lookupOrConstructVector(mesh, "Jre"),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    deltaJim_
    (
        lookupOrConstructVector
        (
            mesh,
            "deltaJim",
            lookupOrConstructVector(mesh, "Jim"),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::harmonicElectromagneticModel> Foam::harmonicElectromagneticModel::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return electromagneticModel::New<harmonicElectromagneticModel>(mesh, phaseName);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::harmonicElectromagneticModel::~harmonicElectromagneticModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::harmonicElectromagneticModel::PotE(bool imaginary)
{
    return imaginary ? PotEim_ : PotEre_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::J(bool imaginary)
{
    return imaginary ? Jim_ : Jre_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::B(bool imaginary)
{
    return imaginary ? Bim_ : Bre_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::deltaJ(bool imaginary)
{
    return imaginary ? deltaJim_ : deltaJre_;
}

//const-access

const Foam::volScalarField& Foam::harmonicElectromagneticModel::PotE(bool imaginary) const
{
    return imaginary ? PotEim_ : PotEre_;
}

const Foam::volVectorField& Foam::harmonicElectromagneticModel::J(bool imaginary) const
{
    return imaginary ? Jim_ : Jre_;
}

const Foam::volVectorField& Foam::harmonicElectromagneticModel::B(bool imaginary) const
{
    return imaginary ? Bim_ : Bre_;
}

const Foam::volVectorField& Foam::harmonicElectromagneticModel::deltaJ(bool imaginary) const
{
    return imaginary ? deltaJim_ : deltaJre_;
}

bool Foam::harmonicElectromagneticModel::isComplex() const
{
    return isComplex_;
}


// ************************************************************************* //
