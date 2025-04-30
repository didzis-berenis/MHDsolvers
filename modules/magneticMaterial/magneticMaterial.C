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

#include "magneticMaterial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(magneticMaterial, 0);
    addToRunTimeSelectionTable(solver, magneticMaterial, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::magneticMaterial::magneticMaterial
(
    fvMesh& mesh
)
:
    solver(mesh),
    electroBase(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::magneticMaterial::~magneticMaterial()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalar Foam::solvers::magneticMaterial::maxDeltaT() const
{
    scalar deltaT = fvModels().maxDeltaT();
    return deltaT;
}


void Foam::solvers::magneticMaterial::preSolve()
{}


void Foam::solvers::magneticMaterial::moveMesh()
{}


void Foam::solvers::magneticMaterial::prePredictor()
{}


void Foam::solvers::magneticMaterial::momentumPredictor()
{}


void Foam::solvers::magneticMaterial::thermophysicalPredictor()
{}


void Foam::solvers::magneticMaterial::pressureCorrector()
{}


void Foam::solvers::magneticMaterial::postCorrector()
{}

void Foam::solvers::magneticMaterial::postSolve()
{}

void Foam::solvers::magneticMaterial::solveElectromagnetics()
{}

// ************************************************************************* //
