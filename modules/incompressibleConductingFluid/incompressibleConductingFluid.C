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

    electroPtr
    (
        electromagneticModel::New(mesh)
    ),

    electro(electroPtr),

    U_old_
    (
        IOobject
        (
            "U_old",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_
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

    rho(rho_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingFluid::~incompressibleConductingFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
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
Foam::tmp<Foam::volVectorField>& Foam::solvers::incompressibleConductingFluid::JUB(bool imaginary)
{
    if (imaginary)
    {
        return JUBim_;
    }
    return JUBre_;
}
*/


void Foam::solvers::incompressibleConductingFluid::solveElectromagnetics();
{
    if (electro.correctElectromagnetics())
    {
        // Update deltaU
        electroPtr->updateDeltaU(U_ - U_old_);
        //Correct current density
        electroPtr->correct();
        //Store old velocity for next update
        U_old_ = U_;
    }
}

void Foam::solvers::incompressibleConductingFluid::setCorrectElectromagnetics()
{
    electroPtr->setCorrectElectromagnetics();
}

// ************************************************************************* //
