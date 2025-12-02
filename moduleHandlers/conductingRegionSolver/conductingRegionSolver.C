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
    runTime_(runTime),
    restartInterval_(runTime.controlDict().lookupOrDefault("restartInterval",0)),
    disableCleanup_(runTime.controlDict().lookupOrDefault("disableCleanup",false)),
    waitInterval(runTime.controlDict().lookupOrDefault("waitInterval",0))
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

    if (isIncompressibleConductingVoF())
    {
        phaseNames_ = Pair<word>(
            getIncompressibleConductingVoFPtr_()->getPhase1Name(),
            getIncompressibleConductingVoFPtr_()->getPhase2Name()
        );

        IOdictionary physicalPropertiesPhase1
        (
            IOobject
            (
                "physicalProperties."+phaseNames_.first(),
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        IOdictionary physicalPropertiesPhase2
        (
            IOobject
            (
                "physicalProperties."+phaseNames_.second(),
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        phaseConductivities_ = Pair<scalar>(
            physicalPropertiesPhase1.lookup<scalar>("sigma"),
            physicalPropertiesPhase2.lookup<scalar>("sigma")
        );

        updateMixtureConductivity();

    }
    else
    {
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
        //get region characteristic size
        characteristicSize_ =
        (
            physicalProperties.found("Lchar") ?
            dimensionedScalar("Lchar",dimLength,physicalProperties) :
            dimensionedScalar("Lchar",dimLength,0)
        ).value();
    }
    /*else 
    {
        FatalIOErrorInFunction(runTime.constant())
            << "physicalProperties missing from "
            << runTime.constant()
            << exit(FatalIOError);
    }*/
    
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
    if (isFluid() || isSolid() || isIncompressibleFluid() || isIncompressibleConductingVoF())
    {
        return dynamic_cast<Foam::electroBase*>(solverPtr_);
    }
    return nullptr;
}

Foam::solvers::incompressibleConductingFluid* Foam::conductingRegionSolver::getIncompressibleFluidPtr_()
{
    if (isIncompressibleFluid())
    {
        return dynamic_cast<Foam::solvers::incompressibleConductingFluid*>(solverPtr_);
    }
    return nullptr;
}

Foam::solvers::conductingFluid* Foam::conductingRegionSolver::getFluidPtr_()
{
    if (isFluid())
    {
        return dynamic_cast<Foam::solvers::conductingFluid*>(solverPtr_);
    }
    return nullptr;
}

Foam::solvers::conductingSolid* Foam::conductingRegionSolver::getSolidPtr_()
{
    if (isSolid())
    {
        return dynamic_cast<Foam::solvers::conductingSolid*>(solverPtr_);
    }
    return nullptr;
}

Foam::solvers::incompressibleConductingVoF* Foam::conductingRegionSolver::getIncompressibleConductingVoFPtr_()
{
    if (isIncompressibleConductingVoF())
    {
        return dynamic_cast<Foam::solvers::incompressibleConductingVoF*>(solverPtr_);
    }
    return nullptr;
}

// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField Foam::conductingRegionSolver::getUpdatedConductivity(bool oldAlpha)
{
    if (isIncompressibleConductingVoF())
    {
        return
        dimensionedScalar
        (
            "sigmaDimensions",
            pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
            1.0
        )*(
            phaseConductivities_.first()*getIncompressibleConductingVoFPtr_()->getAlpha1()
            + phaseConductivities_.second()*getIncompressibleConductingVoFPtr_()->getAlpha2()
        );
    }
    else if (isIncompressibleConductingVoF() && oldAlpha)
    {

        return
        dimensionedScalar
        (
            "sigmaDimensions",
            pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
            1.0
        )*(
            phaseConductivities_.first()*getAlpha1Old()
            + phaseConductivities_.second()*(1-getAlpha1Old())
        );
    }
    else
    {
        return getElectro().sigma();
    }
}

void Foam::conductingRegionSolver::updateMixtureConductivity()
{
    //Pout << "Updating mixture conductivity in region " << name_ << endl;
    if (isIncompressibleConductingVoF())
    {
        volScalarField new_sigma = getUpdatedConductivity();
        //scalar sigmaCutoff = gMax(new_sigma())*SMALL;
        volScalarField& sigma = getElectroBasePtr_()->getSigma();
        forAll(new_sigma, cellI)
        {
            sigma[cellI] = new_sigma[cellI];// < sigmaCutoff ? 0 : new_sigma[cellI];
        }
        //sigma.correctBoundaryConditions();
    }
}

void Foam::conductingRegionSolver::solveElectromagnetics()
{
    updateMixtureConductivity();

    if (getElectroBasePtr_())
    {
        getElectroBasePtr_()->solveElectromagnetics();
    }
}

void Foam::conductingRegionSolver::electromagneticPredictor()
{
    if (getElectroBasePtr_())
    {
        getElectroBasePtr_()->electromagneticPredictor();
    }
}

// Get velocity
const Foam::volScalarField& Foam::conductingRegionSolver::getAlpha1()
{
    if (isIncompressibleConductingVoF())
    {
        return getIncompressibleConductingVoFPtr_()->getAlpha1();
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << incompressibleConductingVoFSolverName_
        << "!\n" << "Cannot get Alpha1 field!\n"
        << exit(FatalIOError);
    }
}

// Get old velocity
const Foam::volScalarField& Foam::conductingRegionSolver::getAlpha1Old()
{
    if (isIncompressibleConductingVoF())
    {
        return getIncompressibleConductingVoFPtr_()->getAlpha1Old();
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << incompressibleConductingVoFSolverName_
        << "!\n" << "Cannot get Alpha1 field!\n"
        << exit(FatalIOError);
    }
}

// Get velocity
const Foam::volVectorField& Foam::conductingRegionSolver::getU()
{
    if (isFluid())
    {
        return getFluid().U;
    }
    else if (isIncompressibleFluid())
    {
        return getIncompressibleFluid().U;
    }
    else if (isIncompressibleConductingVoF())
    {
        return getIncompressibleConductingVoF().U;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << incompressibleFluidSolverName_ << "!\n" << "Cannot get Velocity field!\n"
        << exit(FatalIOError);
    }
}

// Get old velocity
const Foam::volVectorField& Foam::conductingRegionSolver::getUold()
{
    if (isFluid())
    {
        return getFluid().U_old;
    }
    else if (isIncompressibleFluid())
    {
        return getIncompressibleFluid().U_old;
    }
    else if (isIncompressibleConductingVoF())
    {
        return getIncompressibleConductingVoF().U_old;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << incompressibleFluidSolverName_ << "!\n" << "Cannot get Velocity field!\n"
        << exit(FatalIOError);
    }
}

// Get temperature
const Foam::volScalarField& Foam::conductingRegionSolver::getTemperature()
{
    if (isFluid())
    {
        return getFluid().thermo.T();
    }
    else if (isSolid())
    {
        return getSolid().thermo.T();
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Temperature field!\n"
        << exit(FatalIOError);
    }
}

//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolver::setJ(volVectorField& globalField, bool imaginary)
{
    if (getElectroBasePtr_())
    {
        volVectorField& regionField = getElectroBasePtr_()->getJ(imaginary);
        forAll(regionField, cellI)
        {
            regionField[cellI] = globalField[cellI];
        }
        regionField.correctBoundaryConditions();
    }
}
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolver::setB(volVectorField& globalField, bool imaginary)
{
    if (getElectroBasePtr_())
    {
        volVectorField& regionField = getElectroBasePtr_()->getB(imaginary);
        forAll(regionField, cellI)
        {
            regionField[cellI] = globalField[cellI];
        }
        regionField.correctBoundaryConditions();
    }
}
//Assigns U_old = U
void Foam::conductingRegionSolver::storeU()
{
    if (getFluidPtr_())
    {
        getFluidPtr_()->storeU();
    }
    else if(getIncompressibleFluidPtr_())
    {
        getIncompressibleFluidPtr_()->storeU();
    }
}
void Foam::conductingRegionSolver::updateNu(volScalarField& globalField)
{
    if(isIncompressibleConductingVoF())
    {
        globalField = getIncompressibleConductingVoFPtr_()->getNu();
    }
}
void Foam::conductingRegionSolver::storeAlpha1()
{
    if(getIncompressibleConductingVoFPtr_())
    {
        getIncompressibleConductingVoFPtr_()->storeAlpha1();
    }
}
//returns read-only access to electro module
const Foam::electromagneticModel& Foam::conductingRegionSolver::getElectro()
{
    if (getElectroBasePtr_())
    {
        return getElectroBasePtr_()->electro;
    }
    else
    {
        FatalIOError << " electroBase class not found for region "
        << name_ << "!\n" << "Cannot get electromagnetic model!\n" << exit(FatalIOError);
    }
}
const Foam::solvers::incompressibleConductingFluid& Foam::conductingRegionSolver::getIncompressibleFluid()
{
    if (getIncompressibleFluidPtr_())
    {
        const Foam::solvers::incompressibleConductingFluid& fluidRef(*getIncompressibleFluidPtr_());
        return fluidRef;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not "
        << incompressibleFluidSolverName_ << "!\n" << "Cannot get reference!\n"
        << exit(FatalIOError);
    }
}
const Foam::solvers::conductingFluid& Foam::conductingRegionSolver::getFluid()
{
    if (getFluidPtr_())
    {
        const Foam::solvers::conductingFluid& fluidRef(*getFluidPtr_());
        return fluidRef;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not " << fluidSolverName_
        << "!\n" << "Cannot get reference!\n"
        << exit(FatalIOError);
    }
}
const Foam::solvers::conductingSolid& Foam::conductingRegionSolver::getSolid()
{
    if (getSolidPtr_())
    {
        const Foam::solvers::conductingSolid& solidRef(*getSolidPtr_());
        return solidRef;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not "
        << solidSolverName_ << "!\n" << "Cannot get reference!\n"
        << exit(FatalIOError);
    }
}

const Foam::solvers::incompressibleConductingVoF& Foam::conductingRegionSolver::getIncompressibleConductingVoF()
{
    if (getIncompressibleConductingVoFPtr_())
    {
        const Foam::solvers::incompressibleConductingVoF& vofRef(*getIncompressibleConductingVoFPtr_());
        return vofRef;
    }
    else
    {
        FatalIOError
        << " region " << name_ << " solver is not "
        << incompressibleConductingVoFSolverName_ << "!\n" << "Cannot get reference!\n"
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

bool Foam::conductingRegionSolver::isIncompressibleConductingVoF()
{
    return name_ == incompressibleConductingVoFSolverName_;
}

bool Foam::conductingRegionSolver::isElectroHarmonic()
{
    if (isFluid() || isSolid() || isIncompressibleFluid() || isIncompressibleConductingVoF())
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

        if (isFluid() || isIncompressibleFluid())
        {
            const volVectorField& U = getU();
            const volVectorField& U_old = getUold();
            scalar maxMagneticReynoldsDifference = mu_0 * characteristicSize_ *
                gMax((getElectro().sigmaInv()*mag(U_old-U))());
            scalar maxRelativeVelocityDifference =gMax((mag(U_old-U)/(gAverage((mag(U))())+smallU.value()))());
            maxRemDiff_local = max(
                maxMagneticReynoldsDifference,
                maxRemDiff_local);        

            maxRelDiff_local = max(
                maxRelativeVelocityDifference,
                maxRelDiff_local);
        }

        if (isIncompressibleConductingVoF())
        {
            /*volScalarField sigmaNormalized = getElectro().sigma()/
            dimensionedScalar
            (
                "sigmaMag",
                pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
                gMax(getElectro().sigma()())
            );
            volScalarField sigmaOld = getUpdatedConductivity(true);
            volScalarField sigmaNormalizedOld = sigmaOld/
            dimensionedScalar
            (
                "sigmaMagOld",
                pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
                gMax(sigmaOld())
            );
            const volVectorField U = getU()*sigmaNormalized;
            const volVectorField U_old = getUold()*sigmaNormalizedOld;

            scalar maxRelativeVelocityDifference =gMax((mag(U_old-U)/(gAverage((mag(U))())+smallU.value()))());

            maxRelDiff_local = max(
                maxRelativeVelocityDifference,
                maxRelDiff_local);*/
                
            const volScalarField& alpha1 = getAlpha1();
            const volScalarField& alpha1_old = getAlpha1Old();
            scalar maxRelativePhaseDifference =gMax(mag(alpha1_old-alpha1)());

            maxRelDiff_local = max(
                maxRelativePhaseDifference,
                maxRelDiff_local);
            //Info << "maxRelDiff_local: " << maxRelDiff_local
            //<< endl;
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

bool Foam::conductingRegionSolver::needsCleanup()
{
    if (disableCleanup_)
    {
        return false;
    }
    bool doCleanup = true;
    if (restartInterval_ > 0 && cleanupCounter_ % restartInterval_ == 0)
    {
        if (runTime_.startTime().value() > 0 || cleanupCounter_ > 0 )
        {
            // Do not clean initial time step when calculation restarted.
            doCleanup = false;
        }
    }
    return doCleanup;
}

void Foam::conductingRegionSolver::countToCleanup()
{
    cleanupCounter_++;
}

void Foam::conductingRegionSolver::calcTemperatureGradient()
{
    
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime_.constant(),
            solverPtr_->mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    if
    ( 
        physicalProperties.found("temperature_multiplier") &&
        physicalProperties.found("temperature_addition")
    )
    {
        dimensionedVector temperature_multiplier
        (
            "temperature_multiplier",
            dimTemperature/dimLength,
            physicalProperties
        );
        dimensionedScalar temperature_addition
        (
            "temperature_addition",
            dimTemperature,
            physicalProperties
        );

        if (getSolidPtr_())
        {
            volScalarField& T = getSolidPtr_()->getTemperature();
            const volVectorField& coords = T.mesh().C();
            T = (coords & temperature_multiplier) +  temperature_addition;
            /*
            // If more advanced function is necessary, where
            // coordinate dependent function is applied,
            // then iterate through cells.
            forAll (T, cellI)
            {
                vector point = coords[cellI];//x,y,z at cellI
                T[cellI] = temperature_function(coords).value();
            }
            */
            T.write();
        }
        else if (getFluidPtr_())
        {
            volScalarField& T = getFluidPtr_()->getTemperature();
            const volVectorField& coords = T.mesh().C();
            T = (coords & temperature_multiplier) +  temperature_addition;
            T.write();
        }
        /*
        // Display warning message if solver is not fluid or solid
        else
        {
            Info << "Warning: region " << name_ << " solver is not " << fluidSolverName_
            << " or " << solidSolverName_ << "!\n" << "Cannot set temperature gradient!\n"; 
        }
        */
    }
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
