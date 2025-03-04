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

#include "conductingRegionSolvers.H"
#include "Time.H"
#include "globalMeshNew.H"
//#include <experimental/filesystem>
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::conductingRegionSolvers(const Time& runTime)
:
    runTime_(runTime),
    restartInterval_(runTime.controlDict().lookupOrDefault("restartInterval",0)),
    waitInterval(runTime.controlDict().lookupOrDefault("waitInterval",0))
{
    if (runTime.controlDict().found("regionSolvers"))
    {
        const dictionary& conductingRegionSolversDict =
            runTime.controlDict().subDict("regionSolvers");

        forAllConstIter(dictionary, conductingRegionSolversDict, iter)
        {
            const word regionName(iter().keyword());
            const word solverName(iter().stream());

            names_.append(Pair<word>(regionName, solverName));
        }
    }
    else
    {
        FatalIOErrorInFunction(runTime.controlDict())
            << "regionSolvers list missing from "
            << runTime.controlDict().name()
            << exit(FatalIOError);
    }

    regions_.setSize(names_.size());
    solvers_.setSize(names_.size());
    prefixes_.setSize(names_.size());

    string::size_type nRegionNameChars = 0;

    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        const word& solverName = names_[i].second();
        regionIdx_[regionName] = i;

        // Load the solver library
        solver::load(solverName);

        regions_.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        //Foam::autoPtr<Foam::solver> solverPtr = solver::New(solverName, regions_[i]);
        solvers_.set(i, solver::New(solverName, regions_[i]));

        prefixes_[i] = regionName;
        nRegionNameChars = max(nRegionNameChars, regionName.size());
    }
    //Initialize coupled boundary conditions after all solvers have been loaded,
    //because coupledElectricPotentialFvPatchScalarField constructor requires
    //also the solver of the neighbour patch to be loaded.
    //Initialize PotE before deltaJ BCs, because deltaJ BCs require PotE BCs.
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        bool imaginary = isElectroHarmonic();
        evaluatePotEBfs_(regionName);
        if (imaginary)
        {
            evaluatePotEBfs_(regionName,true);
        }
    }
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        bool imaginary = isElectroHarmonic();
        evaluateDeltaJBfs_(regionName);
        if (imaginary)
        {
            evaluateDeltaJBfs_(regionName,true);
        }
        evaluateJBfs_(regionName);
        if (imaginary)
        {
            evaluateJBfs_(regionName,true);
        }
    }

    nRegionNameChars++;

    prefix0_.append(nRegionNameChars, ' ');

    forAll(names_, i)
    {
        prefixes_[i].append(nRegionNameChars - prefixes_[i].size(), ' ');
    }

    //get region characteristic sizes
    characteristicSizes_.setSize(names_.size());
    forAll(names_, i)
    {
        IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                regions_[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        characteristicSizes_[i] =
        (
            physicalProperties.found("Lchar") ?
            dimensionedScalar("Lchar",dimLength,physicalProperties) :
            dimensionedScalar("Lchar",dimLength,0)
        ).value();
    }
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
    // Get local to global cell id mapping
    int globalCellI = 0;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const cellList cells = regionMesh.cells();
        forAll(cells, cellI)
        {
            regionToGlobalCellId[std::make_pair(regionName,cellI)] = globalCellI++;
        }
    }
    // Find regions, which need feedback control
    // And initialize feedback loop controller
    setUpFeedbackControllers_();

    if (!isElectroHarmonic())
    {
        forAll(names_, i)
        {
            const word& regionName = names_[i].first();
            // Get phase shift for all electromagnetic sources
            if (getElectroBasePtr_(regionName))
            {
                electroBase* electrPtr = getElectroBasePtr_(regionName);
                word regionRole = electrPtr->electro.getRegionRole();
                if (regionRole == "coil")
                    conductorPhases_[regionName] = electrPtr->electro.lookup<scalar>("phase");
            }
        }
    }
    checkIfAnyElectricSources_();
    if (hasElectricSources_ && !isElectroHarmonic())
    {
        // Check if frequency is non-negative
        if (frequency_ < 0 )
        {
            FatalIOError << "Frequency must be non-negative, but "
            << frequency_ << "Hz was provided!\n" << exit(FatalIOError);
        }
        // Check if all electromagnetic sources has phase defined
        // for a transient electromagnetic simulation
        forAll(names_, i)
        {
            const word& regionName = names_[i].first();
            if (getElectroBasePtr_(regionName) && getElectroBasePtr_(regionName)->electro.getRegionRole() == "coil")
            {
                bool sourceHasPhase = false;
                for (auto element : conductorPhases_)
                {
                    if (regionName == element.first)
                    {
                        sourceHasPhase = true;
                        break;
                    }
                }
                if (!sourceHasPhase)
                {
                    FatalIOError << "Failed to get phase for electromagnetic source region " << regionName << "!\n"
                    << exit(FatalIOError);
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::~conductingRegionSolvers()
{}

void Foam::conductingRegionSolvers::setUpFeedbackControllers_()
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            electroBase* electrPtr = getElectroBasePtr_(regionName);
            const word feedbackType = electrPtr->electro.lookupOrDefault<word>("feedbackControl","");
            if (feedbackType == "")
                continue;
            if (feedbackType == "current" || feedbackType == "voltage")
            {
                const word terminalName = electrPtr->electro.lookup<word>("terminalName");
                Info << "Setting up controls for terminal: " << terminalName << endl;
                if (terminalToRegions_.find(terminalName) != terminalToRegions_.end())
                {
                    terminalToRegions_[terminalName].push_back(regionName);
                    continue;
                }
                terminalToRegions_[terminalName].push_back(regionName);
                const scalar target_phase =
                dimensionedScalar(
                    "phaseShift",
                    dimless,
                    electrPtr->electro
                ).value();
                if (feedbackType == "current")
                {
                    const scalar target_value =
                    dimensionedScalar(
                        "currentDensity",
                        dimCurrent/dimLength/dimLength,
                        electrPtr->electro
                    ).value();
                    const scalar initial_value = readControlValue_("coilVoltages/"+terminalName);
                    //Info << "Value from file: " << initial_value << endl;
                    //Info << "Value from dict: " << setVoltage.value() << endl;
                    const scalar initial_phase = readControlValue_("coilPhases/"+terminalName);
                    //Info << "Phase from file: " << initial_phase << endl;
                    //Info << "Phase from dict: " << electrPtr->electro.lookup<scalar>("voltagePhase") << endl;
                    /*
                    TODO: Check if the initially guessed proportionality coefficient of target_value isn't grossly incorrect when
                    ferromagnetic materials are present.
                    This could lead to so inefficient (slow) current correction that it doesn't actually work for practical purposes. 
                    This may be avoided with the use of integral regulation.
                    */
                    feedbackControllers_[terminalName] = 
                        feedbackLoopController
                        (
                            Pair<scalar>(target_value,target_phase),
                            feedbackType,
                            Pair<scalar>(initial_value,initial_phase)
                        );
                    scalar phaseReference = 180.0;//PI
                    feedbackControllers_[terminalName].setReference(Pair<scalar>(target_value,phaseReference));
                    // Variables are linked, so phase error is limited by value error
                    scalar maxCurrentError = MAX_CONTROL_ERROR;
                    // Assuming eeror for case Imag=(1-maxCurrentError)*Ire, we can arrive at approximate max error for phase,
                    // where coefficient 1.5 is added, because it is unlikely to reach the theoretical minimum error.
                    scalar maxPhaseError = 1.5*(1.0-maxCurrentError)/sqrt(2.0/maxCurrentError-1.0)/PI;// Divided by PI, because in radians.
                    if (!isElectroHarmonic())
                    {
                        // Additional error is introduced in transient simulations,
                        // whicjh is due to finite time step size.
                        tStepElectro_ = adjustableRunTime_ ?
                        readScalar(runTime_.controlDict().lookup("writeInterval")) :
                        readScalar(runTime_.controlDict().lookup("deltaT"))*writeMultiplier_;
                        frequency_ = runTime_.controlDict().lookup<scalar>("frequency");
                        integrationSteps_ = round(0.5/(tStepElectro_*frequency_));
                        scalar half_step_weight = 1.0/(2.0*integrationSteps_);//Error related to half of step size
                        scalar maxIntegrationError = asin(sin(PI*(1.0-half_step_weight)))/PI;
                        maxPhaseError += maxIntegrationError;
                    }
                    if (maxPhaseError > 0.05)
                    {
                        Info << "Warning: Phase error " << maxPhaseError*100
                        << "% is larger than 5% ! Consider decreasing time step." << endl;
                    }
                    feedbackControllers_[terminalName].setMaxError(Pair<scalar>(maxCurrentError,maxPhaseError));
                    // Stabilizer detects positive feedback loop and stops controller
                    feedbackControllers_[terminalName].setStabilizer(Pair<bool>(false,true),Pair<int>(0,2));
                    feedbackControllers_[terminalName].setMinMaxValue(Pair<scalar>(-great,-phaseReference),Pair<scalar>(great,phaseReference));
                    feedbackControllers_[terminalName].updateCoefficients(
                        Pair<scalar>(
                            target_value != 0 ? -0.5*initial_value/target_value : 0 , 
                            0.5
                        ),
                        Pair<scalar>(0,0),
                        Pair<scalar>(0,0));
                        //Info << "voltage: " << control_value << " " << control_phase <<endl;
                        //Info << "controls: " << control_value/target_value << " " << control_phase - target_phase <<endl;
                }
                if (feedbackType == "voltage")
                {
                    //setVoltage = target_value + induced_voltage
                    //induced_voltage = -0.5*windingDirection*(
                    //innerArea*integralCore(dB/dt)+outerArea*(integralCore(dB/dt)+integralCoil(dB/dt))
                    //) = -0.5*windingDirection*(
                    //(innerArea+outerArea)*integralCore(dB/dt)+outerArea*(integralCoil(dB/dt))
                    //)
                    const scalar target_value =
                    dimensionedScalar(
                        "voltage",
                        dimMass*dimLength*dimLength/(dimTime*dimTime*dimTime*dimCurrent),
                        electrPtr->electro
                    ).value();
                    const scalar inner_area =
                    dimensionedScalar(
                        "innerArea",
                        dimLength*dimLength,
                        electrPtr->electro
                    ).value();
                    const scalar outer_area =
                    dimensionedScalar(
                        "outerArea",
                        dimLength*dimLength,
                        electrPtr->electro
                    ).value();
                    const Foam::vector winding_direction =
                    dimensionedVector(
                        "windingDirection",
                        dimless,
                        electrPtr->electro
                    ).value();
                    coilParameters_ terminalControls;
                    terminalControls.target_value = target_value;
                    terminalControls.inner_area = inner_area;
                    terminalControls.outer_area = outer_area;
                    terminalControls.winding_direction = winding_direction;
                    voltageControls_[terminalName] = terminalControls;
                    //TODO: Get this from terminal boundary
                    scalar set_voltage = target_value;
                    inducedVoltage_[terminalName] = set_voltage - target_value;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
//Assigns single fluid/solid region values from global field to region
void Foam::conductingRegionSolvers::scalarGlobalToField_(volScalarField& global,volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[regionToGlobalCellId[std::make_pair(regionName,cellI)]];
    }
}

void Foam::conductingRegionSolvers::vectorGlobalToField_(volVectorField& global,volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[regionToGlobalCellId[std::make_pair(regionName,cellI)]];
    }
}

Foam::electroBase* Foam::conductingRegionSolvers::getElectroBasePtr_(const word regionName)
{
    if (isFluid(regionName) || isSolid(regionName) || isElectric(regionName))
    {
    Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
        return dynamic_cast<Foam::electroBase*>(basePtr);
    }
    return nullptr;
}

//Note: Could make new solver base class with methods for conducting media
//to avoid casting to each of derived class for accessing new common methods.
Foam::solvers::conductingFluid* Foam::conductingRegionSolvers::getFluidPtr_(const word regionName)
{
    if (isFluid(regionName))
    {
        Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
        return dynamic_cast<Foam::solvers::conductingFluid*>(basePtr);
    }
    return nullptr;
}

Foam::solvers::conductingSolid* Foam::conductingRegionSolvers::getSolidPtr_(const word regionName)
{
    if (isSolid(regionName))
    {
        Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
        return dynamic_cast<Foam::solvers::conductingSolid*>(basePtr);
    }
    return nullptr;
}
//initializes boundary conditions
void Foam::conductingRegionSolvers::evaluatePotEBfs_(const word regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName))
    {
        getElectroBasePtr_(regionName)->initPotE(imaginary);
    }
}
void Foam::conductingRegionSolvers::updatePotErefGrad_(const word regionName, scalar newGrad, bool imaginary)
{
    if (getElectroBasePtr_(regionName))
    {
        getElectroBasePtr_(regionName)->updatePotErefGrad(newGrad, imaginary);
    }
}
//initializes boundary conditions
void Foam::conductingRegionSolvers::evaluateDeltaJBfs_(const word regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName))
    {
        getElectroBasePtr_(regionName)->initDeltaJ(imaginary);
    }
}
void Foam::conductingRegionSolvers::evaluateJBfs_(const word regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName))
    {
        electroBase* electrPtr = getElectroBasePtr_(regionName);
        if (electrPtr->electro.getRegionRole() == "coil")
            electrPtr->initJ(imaginary);
    }
}

void Foam::conductingRegionSolvers::checkIfAnyElectricSources_()
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            const word regionRole = getElectroBasePtr_(regionName)->electro.getRegionRole();
            if (regionRole == "coil")
            {
                // For these regions current density is calculated on OpenFOAM side
                // and incorporated in Elmer as an external current source.
                hasElectricSources_ = true;
                break;
            }
        }
    }
}

Foam::volVectorField Foam::conductingRegionSolvers::getJdirection_(word regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName) && getElectroBasePtr_(regionName)->electro.getRegionRole() == "wire")
    {
        volVectorField regionField = getElectroBasePtr_(regionName)->getJ(imaginary);
        volScalarField magnitudeField = mag(regionField);
        volVectorField directionField = regionField/magnitudeField;
        return directionField;
    }
    else
    {
        FatalIOError << "Failed to get J from region " << regionName << "!\n"
        << exit(FatalIOError);
    }
}

Foam::scalar Foam::conductingRegionSolvers::getCurrentSum_(volVectorField JGlobal,word regionName,bool imaginary)
{
    const volVectorField Jdirection = getJdirection_(regionName,imaginary);//Get reference
    volVectorField JRegion = getElectro(regionName).J(imaginary);//Use copy of internal J for rewriting (has respective regions' mesh)
    vectorGlobalToField_(JGlobal,JRegion,regionName);//Assign received field
    const scalarField Jprojection(JRegion & Jdirection);//Project to a reference
    const scalarField volume = mesh(regionName).V();
    //gSum() (also gAverage()) is multi-threading-safe way to sum over all grid points.
    scalar sumJ = gSum(volume*Jprojection)/max(gSum(volume), vSmall);//Find volume average
    return sumJ;
}

Foam::scalar Foam::conductingRegionSolvers::getTransientPhaseShift_(scalar argument)
{
    scalar angle = 0;
    // Safety check to avoid error, if argument is out of bounds
    if (argument > 1.0)
        angle = 0.5*PI;
    else if (argument < -1.0)
        angle = -0.5*PI;
    else
        angle = asin(argument);
    // Taking midpoint of step gives slightly better performance,
    // because we get the average value over time step, which is
    // not the same as instantaneous value at the present time.
    scalar time_position = runTime_.userTimeValue() - 0.5*tStepElectro_;
    bool wave_condition = static_cast <int> (floor(2*time_position*frequency_)) % 2 == 1;
    bool half_wave_condition = (
        (static_cast <int> (floor(4*time_position*frequency_)) % 2 == 1)
        &&
        (4*time_position*frequency_ > floor(4*time_position*frequency_))
    );
    int phase_sign = (wave_condition ? 1 : -1) * (half_wave_condition ? -1 : 1);
    scalar offset = half_wave_condition ? -PI : 0;
    // Voltage is assigned as voltage = magnitude*Sin(omega*t-phase)
    // Passing argument as voltage/magnitude gives solution for phase as follows.
    scalar time_shift = 2*frequency_*time_position*PI - floor(2*frequency_*time_position)*PI;
    scalar phase_shift = -phase_sign*angle + time_shift + offset;
    return phase_shift;
}

void Foam::conductingRegionSolvers::writeControlValue_(const word fileName, const scalar outputValue)
{
    if (Pstream::master())
    {
        std::ofstream outputFile(fileName, std::ios::out);
        if (outputFile.is_open())
        {
            outputFile << outputValue << std::endl;
            outputFile.close();
        }
        else FatalErrorInFunction << "ERROR: Couldn't open " << fileName << " for writing!\n" << abort(FatalError);
    }
}

Foam::scalar Foam::conductingRegionSolvers::readControlValue_(const word fileName)
{
    scalar inputValue = 0.0;
    std::ifstream inputFile(fileName, std::ios::in);
    if (inputFile.is_open())
    {
        inputFile >> inputValue;
        inputFile.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << fileName << " for reading!\n" << abort(FatalError);
    return inputValue;
}

// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

void Foam::conductingRegionSolvers::setGlobalPrefix() const
{
    Sout.prefix() = prefix0_;
}

void Foam::conductingRegionSolvers::setPrefix(const label i) const
{
    Sout.prefix() = prefixes_[i];
}

void Foam::conductingRegionSolvers::resetPrefix() const
{
    Sout.prefix() = string::null;
}

Foam::fvMesh& Foam::conductingRegionSolvers::mesh(const word regionName)
{
    return regions_[regionIdx_[regionName]];
}

void Foam::conductingRegionSolvers::updateFeedbackControl(volVectorField& JreGlobal,volVectorField& JimGlobal)
{
    if (!isElectroHarmonic())
    {
        integrationCounter_++;
    }
    for (auto element : terminalToRegions_)
    {
        const word terminalName = element.first;
        //Info << "terminalName: " << terminalName;
        if (feedbackControllers_[terminalName].getControlType() == "current")
        {
            scalar avgJre = 0;
            scalar avgJim = 0;
            for (word regionName : element.second)
            {
                //Info << "Region: " << regionName;
                //TODO: check if this works correctly if has multiple regions
                avgJre += getCurrentSum_(JreGlobal,regionName);
                if (isElectroHarmonic())
                {
                    avgJim += getCurrentSum_(JimGlobal,regionName,true);
                }
            }
            
            if (!isElectroHarmonic())
            {
                if (integrationCounter_ % integrationSteps_ != 0)
                {
                    //Info << "Incrementing integralCurrent_ by avgJre " << avgJre << endl;
                    integralCurrent_ += std::abs(avgJre)*tStepElectro_;
                    return;
                }
                integralCurrent_ += std::abs(avgJre)*tStepElectro_;
                integralCurrent_ *= PI*frequency_;
            }
            
            const scalar present_value = isElectroHarmonic() ? std::sqrt(std::pow(avgJre,2)+std::pow(avgJim,2))
            : integralCurrent_;
            const scalar present_phase = isElectroHarmonic() ?
            atan2(avgJim,avgJre)*180/PI :
            getTransientPhaseShift_(avgJre/integralCurrent_)*180/PI;//This is at midpoint of avgJre
            integralCurrent_ = 0;
            Pair<scalar> control_values = 
                feedbackControllers_[terminalName].calculateCorrection
                (
                    Pair<scalar>(present_value,present_phase),
                    runTime_.userTimeValue()
                );
            if (feedbackControllers_[terminalName].needsUpdate())
            {
                Info << "Starting feedback update " << endl;
                Info << "Updating voltage for terminal " << terminalName << endl;
                Info << "Present values: (" << present_value << " " << present_phase << ")" << endl
                << "New control values: " << control_values << endl;
                writeControlValue_("coilVoltages/"+terminalName,control_values.first());//std::abs(control_values.first()));
                writeControlValue_("coilPhases/"+terminalName,control_values.second());
            }
            else
            {
                Info << "Calculated current for terminal " << terminalName  << " is within acceptable errors." << endl;
            }
        }
        if (feedbackControllers_[terminalName].getControlType() == "voltage")
        {
            //TODO: set up correct settings for voltage type
            coilParameters_ terminalControls = voltageControls_[terminalName];
            //TODO: calculate integral magnetic flux change
            //setVoltage = target_value + induced_voltage
            //induced_voltage = -0.5*windingDirection*(
            //innerArea*integralCore(dB/dt)+outerArea*(integralCore(dB/dt)+integralCoil(dB/dt))
            //) = -0.5*windingDirection*(
            //(innerArea+outerArea)*integralCore(dB/dt)+outerArea*(integralCoil(dB/dt))
            //)
            /*terminalControls.target_value = target_value;
            terminalControls.inner_area = inner_area;
            terminalControls.outer_area = outer_area;
            terminalControls.winding_direction = winding_direction;
            inducedVoltage_[terminalName] = ;*/
        }
    }
}

void Foam::conductingRegionSolvers::solveElectromagnetics(const word regionName)
{
    if (getElectroBasePtr_(regionName))
    {
        getElectroBasePtr_(regionName)->solveElectromagnetics();
    }
}

void Foam::conductingRegionSolvers::electromagneticPredictor(const word regionName)
{
    if (getElectroBasePtr_(regionName))
    {
        getElectroBasePtr_(regionName)->electromagneticPredictor();
    }
}
//Returns global mesh
const Foam::fvMesh& Foam::conductingRegionSolvers::globalMesh()
{
    // If pointer empty, create new global mesh.
    if (globalMesh_.empty())
    {
        globalMesh_ = globalMeshNew_();
    }
    // Check if global mesh is already created.
    if (globalMesh_.valid())
    {
        return *globalMesh_;
    }
    else
    {
        FatalIOError << "Failed to get global mesh!\n"
        << exit(FatalIOError);
    }
}
bool Foam::conductingRegionSolvers::hasElectricSources()
{
    return hasElectricSources_;
}
//Get J from electromagnetic source region
Foam::volVectorField Foam::conductingRegionSolvers::getSourceJ(word regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName) && getElectroBasePtr_(regionName)->electro.getRegionRole() == "coil")
    {
        volVectorField regionField = getElectroBasePtr_(regionName)->getJ(imaginary);
        if (!isElectroHarmonic())
        {
            regionField*=Foam::sin(2*PI*frequency_*runTime_.userTimeValue()-conductorPhases_[regionName]);
        }
        return regionField;
    }
    else
    {
        FatalIOError << "Failed to get J from region " << regionName << "!\n"
        << exit(FatalIOError);
    }
}
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolvers::setJ(volVectorField& globalField, bool imaginary)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            volVectorField& regionField = getElectroBasePtr_(regionName)->getJ(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
            // Update boundary values based on boundary conditions.
            regionField.correctBoundaryConditions();
        }
    }
}
//Assigns fluid and solid region values from global field to a specified region field
void Foam::conductingRegionSolvers::setJToRegion(volVectorField& globalField, const word& regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName))
    {
        volVectorField& regionField = getElectroBasePtr_(regionName)->getJ(imaginary);
        vectorGlobalToField_(globalField,regionField,regionName);
        // Update boundary values based on boundary conditions.
        regionField.correctBoundaryConditions();
    }
}
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolvers::setB(volVectorField& globalField, bool imaginary)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            volVectorField& regionField = getElectroBasePtr_(regionName)->getB(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
            // Update boundary values based on boundary conditions.
            regionField.correctBoundaryConditions();
        }
    }
}
//Assigns U_old = U
void Foam::conductingRegionSolvers::storeU(const word regionName)
{
    if (getFluidPtr_(regionName))
    {
        getFluidPtr_(regionName)->storeU();
    }
}
//returns read-only access to electro module
const Foam::electromagneticModel& Foam::conductingRegionSolvers::getElectro(const word regionName)
{
    if (getElectroBasePtr_(regionName))
    {
        return getElectroBasePtr_(regionName)->electro;
    }
    else
    {
        FatalIOError << " electroBase class not found for region "
        << regionName << "!\n" << "Cannot get electromagnetic model!\n" << exit(FatalIOError);
    }
}
const Foam::solvers::conductingFluid& Foam::conductingRegionSolvers::getFluid(const word regionName)
{
    if (getFluidPtr_(regionName))
    {
        const Foam::solvers::conductingFluid& fluidRef(*getFluidPtr_(regionName));
        return fluidRef;
    }
    else
    {
        FatalIOError << " region " << regionName << " solver is not "
        << fluidSolverName_ << "!\n" << "Cannot get fluid!\n" << exit(FatalIOError);
    }
}
const Foam::solvers::conductingSolid& Foam::conductingRegionSolvers::getSolid(const word regionName)
{
    if (getSolidPtr_(regionName))
    {
        const Foam::solvers::conductingSolid& solidRef(*getSolidPtr_(regionName));
        return solidRef;
    }
    else
    {
        FatalIOError << " region " << regionName << " solver is not "
        << solidSolverName_ << "!\n" << "Cannot get solid!\n" << exit(FatalIOError);
    }
}

//Assigns single fluid/solid region values from region to global field
void Foam::conductingRegionSolvers::scalarFieldToGlobal(volScalarField& global,const volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[regionToGlobalCellId[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}

void Foam::conductingRegionSolvers::vectorFieldToGlobal(volVectorField& global,const volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[regionToGlobalCellId[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}
bool Foam::conductingRegionSolvers::isFluid(const word regionName)
{
    return names_[regionIdx_[regionName]].second() == fluidSolverName_;
}

bool Foam::conductingRegionSolvers::isSolid(const word regionName)
{
    return names_[regionIdx_[regionName]].second() == solidSolverName_;
}

bool Foam::conductingRegionSolvers::isElectric(const word regionName)
{
    return names_[regionIdx_[regionName]].second() == electricSolverName_;
}

bool Foam::conductingRegionSolvers::isSource(const word regionName)
{
    const word regionRole = getElectro(regionName).getRegionRole();
    return isElectric(regionName) && (regionRole == "coil");
}

bool Foam::conductingRegionSolvers::isSolvedFor(const word regionName)
{
    const word regionRole = getElectro(regionName).getRegionRole();
    return isElectric(regionName) && (regionRole == "coil" || regionRole == "wire");
}

bool Foam::conductingRegionSolvers::hasAnyRole(const word regionRole)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectro(regionName).getRegionRole() == regionRole)
        {
            return true;
        }
    }
    return false;
}

void Foam::conductingRegionSolvers::setPotentialCorrectors(const dictionary& dict)
{
    nPotEcorr_ = dict.lookupOrDefault<label>("nPotECorrectors", nPotEcorr_);
}

bool Foam::conductingRegionSolvers::correctElectroPotential()
{
    if (PotEcorr_ >= nPotEcorr_)
    {
        PotEcorr_ = 0;
        return false;
    }

    PotEcorr_++;

    return true;
}

bool Foam::conductingRegionSolvers::isElectroHarmonic()
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (isFluid(regionName) || isSolid(regionName) || isElectric(regionName))
        {
            return getElectro(regionName).isComplex();
        }
    }
    return false;
}

bool Foam::conductingRegionSolvers::updateMagneticField()
{
    bool doUpdate = false;
    if (isElectroHarmonic())
    {
        scalar maxRemDiff_local = SMALL;        
        scalar maxRelDiff_local = SMALL;

        forAll(names_, i)
        {
            const word& regionName = names_[i].first();
            if (isFluid(regionName))
            {
                const volVectorField& U = getFluid(regionName).U;
                const volVectorField& U_old = getFluid(regionName).U_old;
                maxRemDiff_local = max(
                    mu_0 * characteristicSizes_[i] *
                    max(getElectro(regionName).sigmaInv()*mag(U_old-U)).value(),
                    maxRemDiff_local);        

                maxRelDiff_local = max(
                    (max(mag(U_old-U)/(average(mag(U))+smallU))).value(),
                    maxRemDiff_local);
            }
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

bool Foam::conductingRegionSolvers::needsCleanup()
{
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

void Foam::conductingRegionSolvers::countToCleanup()
{
    cleanupCounter_++;
}

void Foam::conductingRegionSolvers::calcTemperatureGradient(const word regionName)
{
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime_.constant(),
            mesh(regionName),
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

        if (getSolidPtr_(regionName))
        {
            volScalarField& T = getSolidPtr_(regionName)->getTemperature();
            T = (mesh(regionName).C() & temperature_multiplier) +  temperature_addition;
            T.write();
        }
        else if (getFluidPtr_(regionName))
        {
            volScalarField& T = getFluidPtr_(regionName)->getTemperature();
            T = (mesh(regionName).C() & temperature_multiplier) +  temperature_addition;
            T.write();
        }
        /*
        // Display warning message if solver is not fluid or solid
        else
        {
            Info << "Warning: region " << regionName << " solver is not " << fluidSolverName_
            << " or " << solidSolverName_ << "!\n" << "Cannot set temperature gradient!\n"; 
        }
        */
    }
}

int Foam::conductingRegionSolvers::getRegionId(const word regionName)
{
    return regionIdx_[regionName];
}

Foam::List<Foam::Pair<Foam::word>> Foam::conductingRegionSolvers::getNames()
{
    return names_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::solver& Foam::conductingRegionSolvers::operator[](const label i)
{
    setPrefix(i);
    return solvers_[i];
}

Foam::solver& Foam::conductingRegionSolvers::operator()(const word regionName)
{
    setPrefix(regionIdx_[regionName]);
    return solvers_[regionIdx_[regionName]];
}
/*
Foam::solver* Foam::conductingRegionSolvers::operator()(const label i)
{
    return solvers_(i);
}*/


// ************************************************************************* //
