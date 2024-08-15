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
:
    runTime_(runTime)
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
    
    //get region characteristic size
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    characteristicSize_ = physicalProperties.lookupOrDefault<dimensionedScalar>
    (
        "Lchar",
        dimensionedScalar
        (
            dimLength,
            0
        )
    ).value();
    
    // get write controls / magnetic field update controls
    if (isElectroHarmonic())
    {
        // Maximum allowable magnetic Reynolds number difference comparing
        // to last magnetic field update.
        // This option controls frequency magnetic field updating is called.
        //     (0,inf) - when relative difference in any cell exceeds given value
        //     0     - every iteration
        maxRemDiff_ = readScalar(runTime.controlDict().lookup("maxRemDiff"));

        // Maximum allowable relative field difference in any cell comparing
        // to last magnetic field update.
        // This option controls frequency magnetic field updating is called.
        //     >1  - once
        //     1 - magnetic Reynolds number is used instead
        //     (0,1) - when relative difference in any cell exceeds given value
        //     0     - every iteration
        maxRelDiff_ = readScalar(runTime.controlDict().lookup("maxRelDiff"));
    }
    else
    {
        writeMultiplier_ = readScalar(runTime.controlDict().lookup("writeMultiplier"));
        writeControlDict_ = word(runTime.controlDict().lookup("writeControl"));
        adjustableRunTime_ = (writeControlDict_=="adjustableRunTime");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductingRegionSolver::~conductingRegionSolver()
{}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::electroBase* Foam::conductingRegionSolver::getElectroBasePtr_()
{
    if (isFluid() || isSolid() || isIncompressibleFluid())
    {
        return dynamic_cast<Foam::electroBase*>(solverPtr_);
    }
    return nullptr;
}

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

// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

void Foam::conductingRegionSolver::solveElectromagnetics()
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_();
    if (electroBasePtr)
    {
        electroBasePtr->solveElectromagnetics();
    }
    /**/
    else
    {
        Info << "Warning: region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot solve electromagnetics!\n"; 
    }
    /**/
}

void Foam::conductingRegionSolver::electromagneticPredictor()
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_();
    if (electroBasePtr)
    {
        electroBasePtr->electromagneticPredictor();
    }
    /**/
    else
    {
        Info << "Warning: region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot predict electromagnetics!\n"; 
    }
    /**/
}

/*
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
*/
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolver::setJ(volVectorField& globalField, bool imaginary)
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_();
    if (electroBasePtr)
    {
        volVectorField& regionField = electroBasePtr->getJ(imaginary);
        forAll(regionField, cellI)
        {
            regionField[cellI] = globalField[cellI];
        }
    }
}
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolver::setB(volVectorField& globalField, bool imaginary)
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_();
    if (electroBasePtr)
    {
        volVectorField& regionField = electroBasePtr->getB(imaginary);
        forAll(regionField, cellI)
        {
            regionField[cellI] = globalField[cellI];
        }
    }
}
//Assigns U_old = U
void Foam::conductingRegionSolver::storeU()
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_();
    if (fluidPtr)
    {
        fluidPtr->storeU();
    }
}
//returns read-only access to electro module
const Foam::electromagneticModel& Foam::conductingRegionSolver::getElectro()
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_();
    if (electroBasePtr)
    {
        return electroBasePtr->electro;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Electric module!\n"
        << exit(FatalIOError);
    }
}
const Foam::solvers::incompressibleConductingFluid& Foam::conductingRegionSolver::getIncompressibleFluid()
{
    Foam::solvers::incompressibleConductingFluid* fluidPtr = getIncompressibleFluidPtr_();
    if (fluidPtr)
    {
        const Foam::solvers::incompressibleConductingFluid& fluidRef(*fluidPtr);
        return fluidRef;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Electric module!\n"
        << exit(FatalIOError);
    }
}
const Foam::solvers::conductingFluid& Foam::conductingRegionSolver::getFluid()
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_();
    if (fluidPtr)
    {
        const Foam::solvers::conductingFluid& fluidRef(*fluidPtr);
        return fluidRef;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Electric module!\n"
        << exit(FatalIOError);
    }
}
const Foam::solvers::conductingSolid& Foam::conductingRegionSolver::getSolid()
{
    Foam::solvers::conductingSolid* solidPtr = getSolidPtr_();
    if (solidPtr)
    {
        const Foam::solvers::conductingSolid& solidRef(*solidPtr);
        return solidRef;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Electric module!\n"
        << exit(FatalIOError);
    }
}

bool Foam::conductingRegionSolver::isIncompressibleFluid()
{
    return name_ == incompressibleFluidSolverName_;
}

bool Foam::conductingRegionSolver::isFluid()
{
    return name_ == fluidSolverName_;
}

bool Foam::conductingRegionSolver::isSolid()
{
    return name_ == solidSolverName_;
}

void Foam::conductingRegionSolver::setCorrectElectromagnetics()
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_();
    if (electroBasePtr)
    {
        electroBasePtr->setCorrectElectromagnetics();
    }
}

bool Foam::conductingRegionSolver::isElectroHarmonic()
{
    if (isFluid() || isSolid() || isIncompressibleFluid())
    {
        return getElectro().isComplex();
    }
    return false;
}

bool Foam::conductingRegionSolver::updateMagneticField()
{
    bool doUpdate = false;
    if (isElectroHarmonic())
    {
        scalar maxRemDiff_local = SMALL;        
        scalar maxRelDiff_local = SMALL;

        if (isFluid())
        {
            const volVectorField& U = getFluid().U;
            const volVectorField& U_old = getFluid().U_old;
            maxRemDiff_local = max(
                mu_0 * characteristicSize_ *
                max(getElectro().sigmaInv()*mag(U_old-U)).value(),
                maxRemDiff_local);        

            maxRelDiff_local = max(
                (max(mag(U_old-U)/(average(mag(U))+smallU))).value(),
                maxRemDiff_local);
        }

        if((maxRelDiff_local>maxRelDiff_ || maxRelDiff_<SMALL) && maxRelDiff_+SMALL<=1.0) {
            doUpdate = true;
        }
        else if(maxRemDiff_local>maxRemDiff_ && maxRelDiff_-SMALL<=1.0) {
            doUpdate = true;
        }
    }
    else
    {
        if (adjustableRunTime_)
        {
            if (runTime_.writeTime()) doUpdate = true;
        }
        else
        {
            writeCounter_++;
            if ( (writeCounter_ % writeMultiplier_) == 0 && runTime_.run()) doUpdate = true;
        }
    }
    return doUpdate;
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
