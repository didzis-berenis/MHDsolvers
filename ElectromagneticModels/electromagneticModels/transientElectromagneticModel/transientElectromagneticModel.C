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
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::transientElectromagneticModel::transientElectromagneticModel
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    electromagneticModel(mesh),
    PotE_(lookupOrConstructScalar(mesh, "PotE")),
    Jre_(lookupOrConstructVector(mesh, "J")),
    Bre_(lookupOrConstructVector(mesh, "B")),
    //Get boundary conditions from J
    deltaJre_
    (
        lookupOrConstructVector
        (
            mesh,
            "deltaJ",
            lookupOrConstructVector(mesh, "J")
        )
    ),
    deltaUxBre_
    (
        lookupOrConstructVector
        (
            mesh,
            "deltaUxB",
            lookupOrConstructVector(mesh, "B")*dimensionedScalar(dimVelocity,0)
        )
    )
{}

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
    return Jre_;
}

Foam::volVectorField& Foam::transientElectromagneticModel::B(bool imaginary)
{
    return Bre_;
}

Foam::volVectorField& Foam::transientElectromagneticModel::deltaJ(bool imaginary)
{
    return deltaJre_;
}

//const-access

const Foam::volScalarField& Foam::transientElectromagneticModel::PotE(bool imaginary) const
{
    return PotE_;
}

const Foam::volVectorField& Foam::transientElectromagneticModel::J(bool imaginary) const
{
    return Jre_;
}

const Foam::volVectorField& Foam::transientElectromagneticModel::B(bool imaginary) const
{
    return Bre_;
}

const Foam::volVectorField& Foam::transientElectromagneticModel::deltaUxB(bool imaginary) const
{
    return deltaUxBre_;
}

void Foam::transientElectromagneticModel::updateDeltaU(volVectorField& Udiff)
{
    deltaUxBre_ = Udiff ^ Bre_;
}

const Foam::volVectorField& Foam::transientElectromagneticModel::deltaJ(bool imaginary) const
{
    return deltaJre_;
}

bool Foam::transientElectromagneticModel::isComplex() const
{
    return isComplex_;
}

// ************************************************************************* //
