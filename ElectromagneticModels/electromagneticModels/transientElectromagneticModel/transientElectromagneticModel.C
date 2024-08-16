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

#include "transientElectromagneticModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(transientElectromagneticModel, 0);
    // Add to runtime selection Table
    /*addToRunTimeSelectionTable
    (
        electromagneticModel,
        transientElectromagneticModel,
        fvMesh
    );*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transientElectromagneticModel::transientElectromagneticModel
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    electromagneticModel(mesh),
    PotE_(lookupOrConstructScalar(mesh, "PotE"))
{
    constructVector(mesh, JName_);
    constructVector(mesh, BName_);
    //Get boundary conditions from J
    deltaJ_ = tmp<volVectorField>(J());
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::transientElectromagneticModel> Foam::transientElectromagneticModel::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return electromagneticModel::New<transientElectromagneticModel>(mesh, phaseName);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::transientElectromagneticModel::~transientElectromagneticModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::transientElectromagneticModel::PotE(bool imaginary)
{
    return PotE_;
}

Foam::volVectorField& Foam::transientElectromagneticModel::J(bool imaginary)
{
    return getVectorFieldRef(JName_);
}

Foam::volVectorField& Foam::transientElectromagneticModel::B(bool imaginary)
{
    return getVectorFieldRef(BName_);
}

Foam::tmp<Foam::volVectorField>& Foam::transientElectromagneticModel::deltaJ(bool imaginary)
{
    return deltaJ_;
}

//const-access

const Foam::volScalarField& Foam::transientElectromagneticModel::PotE(bool imaginary) const
{
    return PotE_;
}

const Foam::volVectorField& Foam::transientElectromagneticModel::J(bool imaginary) const
{
    return getVectorField(JName_);
}

const Foam::volVectorField& Foam::transientElectromagneticModel::B(bool imaginary) const
{
    return getVectorField(BName_);
}

const Foam::tmp<Foam::volVectorField>& Foam::transientElectromagneticModel::deltaJ(bool imaginary) const
{
    return deltaJ_;
}

bool Foam::transientElectromagneticModel::isComplex() const
{
    return isComplex_;
}

// ************************************************************************* //
