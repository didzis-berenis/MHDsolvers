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
    PotEim_(lookupOrConstructScalar(mesh, "PotEim"))
{
    constructVector(mesh, JreName_);
    constructVector(mesh, JimName_);
    constructVector(mesh, BreName_);
    constructVector(mesh, BimName_);
    //Get boundary conditions from J
    deltaJ() = tmp<volVectorField>(J());
    deltaJ(true) = tmp<volVectorField>(J(true));
}


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
    if (imaginary)
    {
        return PotEim_;
    }
    return PotEre_;
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::J(bool imaginary)
{
    if (imaginary)
    {
        return getVectorFieldRef(JimName_);
    }
    return getVectorFieldRef(JreName_);
}

Foam::volVectorField& Foam::harmonicElectromagneticModel::B(bool imaginary)
{
    if (imaginary)
    {
        return getVectorFieldRef(BimName_);
    }
    return getVectorFieldRef(BreName_);
}

Foam::tmp<Foam::volVectorField>& Foam::harmonicElectromagneticModel::deltaJ(bool imaginary)
{
    if (imaginary)
    {
        return deltaJim_;
    }
    return deltaJre_;
}

//const-access

const Foam::volScalarField& Foam::harmonicElectromagneticModel::PotE(bool imaginary) const
{
    if (imaginary)
    {
        return PotEim_;
    }
    return PotEre_;
}

const Foam::volVectorField& Foam::harmonicElectromagneticModel::J(bool imaginary) const
{
    if (imaginary)
    {
        return getVectorField(JimName_);
    }
    return getVectorField(JreName_);
}

const Foam::volVectorField& Foam::harmonicElectromagneticModel::B(bool imaginary) const
{
    if (imaginary)
    {
        return getVectorField(BimName_);
    }
    return getVectorField(BreName_);
}

const Foam::tmp<Foam::volVectorField>& Foam::harmonicElectromagneticModel::deltaJ(bool imaginary) const
{
    if (imaginary)
    {
        return deltaJim_;
    }
    return deltaJre_;
}

bool Foam::harmonicElectromagneticModel::isComplex() const
{
    return isComplex_;
}


// ************************************************************************* //
