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

#include "conductingRegionSolver.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductingRegionSolver::conductingRegionSolver(const Time& runTime, fvMesh& mesh)
{

    if (runTime.controlDict().found("solver"))
    {
        // Read the name_ from the optional solver entry in controlDict
        name_ = word
        (
            runTime.controlDict().lookup("solver")
        );
    }
    else
    {
        FatalIOErrorInFunction(runTime.controlDict())
            << "solver missing from "
            << runTime.controlDict().name()
            << exit(FatalIOError);
    }

    // Load the solver library
    solver::load(name_);

    // Instantiate the selected solver
    Foam::autoPtr<solver> solverAutoPtr(solver::New(name_, mesh));
    // get pointer to solver
    solverPtr_ = solverAutoPtr.ptr();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductingRegionSolver::~conductingRegionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingFluid* Foam::conductingRegionSolver::getIncompressibleFluidPtr_()
{
    if (name_ == incompressibleFluidSolverName_)
    {
        return dynamic_cast<Foam::solvers::incompressibleConductingFluid*>(solverPtr_);
    }
    return nullptr;
}

Foam::solvers::conductingFluid* Foam::conductingRegionSolver::getFluidPtr_()
{
    if (name_ == fluidSolverName_)
    {
        return dynamic_cast<Foam::solvers::conductingFluid*>(solverPtr_);
    }
    return nullptr;
}

Foam::solvers::conductingSolid* Foam::conductingRegionSolver::getSolidPtr_()
{
    if (name_ == solidSolverName_)
    {
        return dynamic_cast<Foam::solvers::conductingSolid*>(solverPtr_);
    }
    return nullptr;
}

// Return region fvMesh
void Foam::conductingRegionSolver::setJxB(Foam::volVectorField& field)
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_();
    Foam::solvers::incompressibleConductingFluid* incompressibleFluidPtr = getIncompressibleFluidPtr_();
    if (fluidPtr)
    {
        fluidPtr->setJxB(field);
    }
    else if (incompressibleFluidPtr)
    {
        incompressibleFluidPtr->setJxB(field);
    }
    else
    {
        Info << "Warning: region " << name_ << " solver is not " << fluidSolverName_
        << " or " << incompressibleFluidSolverName_ << "!\n" << "Cannot set JxB field!\n"; 
    }
}

// Return region fvMesh
void Foam::conductingRegionSolver::setJJsigma(Foam::volScalarField& field)
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_();
    Foam::solvers::conductingSolid* solidPtr = getSolidPtr_();
    if (fluidPtr)
    {
        fluidPtr->setJJsigma(field);
    }
    else if (solidPtr)
    {
        solidPtr->setJJsigma(field);
    }
    else
    {
        Info << "Warning: region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot set JJsigma field!\n"; 
    }
}

// Return region fvMesh
Foam::volVectorField& Foam::conductingRegionSolver::getVelocity()
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_();
    Foam::solvers::incompressibleConductingFluid* incompressibleFluidPtr = getIncompressibleFluidPtr_();
    if (fluidPtr)
    {
        return fluidPtr->getVelocity();
    }
    else if (incompressibleFluidPtr)
    {
        return incompressibleFluidPtr->getVelocity();
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << incompressibleFluidSolverName_ << "!\n" << "Cannot get Velocity field!\n"
        << exit(FatalIOError);
    }
}

// Return region fvMesh
Foam::volScalarField& Foam::conductingRegionSolver::getPressure()
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_();
    Foam::solvers::incompressibleConductingFluid* incompressibleFluidPtr = getIncompressibleFluidPtr_();
    if (fluidPtr)
    {
        return fluidPtr->getPressure();
    }
    else if (incompressibleFluidPtr)
    {
        return incompressibleFluidPtr->getPressure();
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << incompressibleFluidSolverName_ << "!\n" << "Cannot get Pressure field!\n"
        << exit(FatalIOError);
    }
}

// Return region fvMesh
Foam::volScalarField& Foam::conductingRegionSolver::getTemperature()
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_();
    Foam::solvers::conductingSolid* solidPtr = getSolidPtr_();
    if (fluidPtr)
    {
        return fluidPtr->getTemperature();
    }
    else if (solidPtr)
    {
        return solidPtr->getTemperature();
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Temperature field!\n"
        << exit(FatalIOError);
    }
}

bool Foam::conductingRegionSolver::isIncompressibleFluid()
{
    if (name_ == incompressibleFluidSolverName_)
    {
        return true;
    }
    return false;
}

bool Foam::conductingRegionSolver::isFluid()
{
    if (name_ == fluidSolverName_)
    {
        return true;
    }
    return false;
}

bool Foam::conductingRegionSolver::isSolid()
{
    if (name_ == solidSolverName_)
    {
        return true;
    }
    return false;
}

Foam::word Foam::conductingRegionSolver::getName()
{
    return name_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
Foam::solver& Foam::conductingRegionSolver::operator()()
{
    return solver_();
}*/


// ************************************************************************* //
