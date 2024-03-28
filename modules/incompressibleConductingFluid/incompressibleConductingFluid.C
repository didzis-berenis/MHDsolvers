/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "incompressibleConductingFluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleConductingFluid, 0);
    addToRunTimeSelectionTable(solver, incompressibleConductingFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingFluid::incompressibleConductingFluid
(
    fvMesh& mesh
)
:
    incompressibleFluid(mesh),

    JxB_
    (
        IOobject
        (
            "JxB",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimensionSet(1, -2, -2, 0, 0, 0, 0), Foam::vector(0,0,0))
    ),
    rho_
    (
        "rho",
        dimDensity,
        IOdictionary
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    ),
    JxB(JxB_),
    rho(rho_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingFluid::~incompressibleConductingFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleConductingFluid::setJxB(volVectorField& JxB)
{   
    JxB_=JxB;
}

Foam::volVectorField& Foam::solvers::incompressibleConductingFluid::getVelocity()
{   
    return U_;
}

Foam::volScalarField& Foam::solvers::incompressibleConductingFluid::getPressure()
{   
    return p_;
}


// ************************************************************************* //
