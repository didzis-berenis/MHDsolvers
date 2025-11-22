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
#include "triPointRef.H"
#include "globalMeshNew.H"
using std::abs;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::conductingRegionSolvers(const Time& runTime)
:
    runTime_(runTime),
    frequency_(runTime_.controlDict().lookup<scalar>("frequency")),
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

        solvers_.set(i, solver::New(solverName, regions_[i]));

        prefixes_[i] = regionName;
        nRegionNameChars = max(nRegionNameChars, regionName.size());
    }
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        Pout << "regionName: " << regionName << endl;
        Pout << "mesh(regionName).cells().size(): " << mesh(regionName).cells().size() << endl;
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
        const word& regionName = names_[i].first();
        if (isFluid(regionName) || isSolid(regionName))
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
        else
        {
            characteristicSizes_[i] = 0;
        }
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
    checkIfAnyVoltageControl_();
    // Find regions, which need feedback control
    // And initialize feedback loop controller
    setUpFeedbackControllers_();
    for (auto element : terminalToRegions_)
    {
        // Coil regions
        const word terminalName = element.first;
        Info << "terminalName" << terminalName << endl;
        for (word regionName : element.second)
        {
            Info << "regionName" << regionName << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::~conductingRegionSolvers()
{}

void Foam::conductingRegionSolvers::setUpFeedbackControllers_()
{
    if (hasVoltageControl_)
    {
        target_value_error_ = 0.01;// 1%
        target_phase_error_ = 0.05;// 5%
        // Relaxation factor may need to be smaller for cases with
        // high-frequency, high-mu ferromagnetic core or large number of windings.
        induced_voltage_relaxation_factor_ = 0.4;
    }
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            electroBase* electrPtr = getElectroBasePtr_(regionName);
            const word feedbackType = electrPtr->electro.lookupOrDefault<word>("feedbackControl","");
            if (feedbackType == "current" || feedbackType == "voltage")
            {
                const word terminalName = electrPtr->electro.lookup<word>("terminalName");
                // Not tested for multiple regions,
                // but also presently not possible configuration to have several regions for one wire.
                terminalToRegions_[terminalName].push_back(regionName);
            }
            if (hasVoltageControl_)
            {
                if (!isElectroHarmonic())
                {
                    const vectorField oldB = getElectro(regionName).B().internalField();
                    regionOldB_[regionName]=oldB;
                }
            }
        }
    }
    if (!isElectroHarmonic())
    {
        tStepElectro_ = adjustableRunTime_ ?
        readScalar(runTime_.controlDict().lookup("writeInterval")) :
        readScalar(runTime_.controlDict().lookup("deltaT"))*writeMultiplier_;
        integrationSteps_ = round(0.5/(tStepElectro_*frequency_));
    }
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
                const scalar target_phase =
                dimensionedScalar(
                    "phaseShift",
                    dimless,
                    electrPtr->electro
                ).value();
                if (feedbackType == "current")
                {
                    integralCurrent_[terminalName] = 0;// Used to store integral values
                    const scalar target_value =
                    dimensionedScalar(
                        "currentDensity",
                        dimCurrent/dimLength/dimLength,
                        electrPtr->electro
                    ).value();
                    const scalar initial_value = readControlValue_("coilVoltages/"+terminalName);
                    const scalar initial_phase = readControlValue_("coilPhases/"+terminalName);
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
                        // which is due to finite time step size.
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
                }
                if (feedbackType == "voltage")
                {
                    electroBase* electrPtr = getElectroBasePtr_(regionName);
                    IOdictionary electroDict = electrPtr->electro;
                    // Wire region boundary mesh
                    const fvBoundaryMesh& boundaryMesh = mesh(regionName).boundary();

                    //Check if initialized
                    if (Pstream::parRun() && !electroDict.found("windingDirection"))
                    {
                        FatalIOError << "Failed to get voltage parameters for region " << regionName
                        << "! Please initialize with a single thread.\n"
                        << exit(FatalIOError);
                    }
                    if (!Pstream::parRun() && !electroDict.found("windingDirection"))//isDict("wireProperties"))
                    {
                        // The following proceedure is much easier, when run with a single thread.
                        // So it is expected to be first initialized on a single thread.
                        // The found wire properties are appended in the dictionary
                        // and can be simply read when using multi-threading.
                        // 1) Find coilCenter position.
                        Foam::vector coilCenter = getCenter(regionName);
                        Pout << "region: " << regionName << " weighted center: " << coilCenter << endl;
                        // 2) Find terminalCenter position.
                        Foam::vector terminalCenter = getCenter(regionName,terminalName);
                        Pout << "terminal: " << terminalName << " weighted center: " << terminalCenter << endl;
                        // 3) Find groundCenter position.
                        word groundTerminalName = electrPtr->findGroundTerminal(terminalName);
                        Pout << "groundTerminalName: " << groundTerminalName << endl;
                        Foam::vector groundCenter = getCenter(regionName,groundTerminalName);
                        Pout << "groundTerminal: " << groundTerminalName << " weighted center: " << groundCenter << endl;
                        // 4) Construct vector terminalDirection from coilCenter to terminalCenter.
                        Foam::vector terminalDirection = terminalCenter - coilCenter;
                        // 5) Construct vector vg from terminalCenter to groundCenter.
                        Foam::vector vg = terminalCenter - groundCenter;
                        // 6) Calculate approximate windingDirectionAP as a cross product between terminalDirection and vg.
                        Foam::vector windingDirectionAP = (terminalDirection ^ vg)/mag(terminalDirection ^ vg);
                        // 7) From min/max distance to coil center determine inner and outer coil boundary
                        label terminalLabel = -1;
                        forAll(boundaryMesh,patchI)
                        {
                            if (boundaryMesh[patchI].name() == terminalName)
                                terminalLabel = patchI;
                        }
                        //Pout << "terminalLabel: " << terminalLabel << endl;
                        const fvPatch& terminalPatch = boundaryMesh[terminalLabel];
                        const pointField terminalSurfacePoints = terminalPatch.patch().localPoints();
                        label minLabel = -1;
                        scalar minDistance = pow(10,10);
                        label maxLabel = -1;
                        scalar maxDistance = pow(10,10);
                        forAll(terminalSurfacePoints,pointI)
                        {
                            point thisPoint = terminalSurfacePoints[pointI];
                            scalar distance = mag(coilCenter - thisPoint);
                            if (distance < minDistance)
                            {
                                minLabel = pointI;
                                minDistance = distance;
                            }
                        }
                        // Easier to find closest point, so using mirror point
                        Foam::vector mirrorCenter = coilCenter+2*terminalDirection;
                        forAll(terminalSurfacePoints,pointI)
                        {
                            point thisPoint = terminalSurfacePoints[pointI];
                            scalar distance = mag(mirrorCenter - thisPoint);
                            if (distance < maxDistance)
                            {
                                maxLabel = pointI;
                                maxDistance = distance;
                            }
                        }
                        /*Pout << "minLabel: " << minLabel << endl;
                        Pout << "minDistance: " << minDistance << endl;
                        Pout << "maxLabel: " << maxLabel << endl;
                        Pout << "maxDistance: " << maxDistance << endl;*/
                        word innerBoundaryName = "";
                        word outerBoundaryName = "";
                        label innerBoundaryLabel = -1;
                        label outerBoundaryLabel = -1;
                        forAll(boundaryMesh,patchI)
                        {
                            const word thisPatchName = boundaryMesh[patchI].name();
                            if (thisPatchName == terminalName || thisPatchName == groundTerminalName)
                                continue;
                            const pointField thisSurfacePoints = boundaryMesh[patchI].patch().localPoints();
                            forAll(thisSurfacePoints,pointI)
                            {
                                // Find common points between patches
                                if(thisSurfacePoints[pointI] == terminalSurfacePoints[minLabel])
                                {
                                    innerBoundaryName = thisPatchName;
                                    innerBoundaryLabel = patchI;
                                }
                                if(thisSurfacePoints[pointI] == terminalSurfacePoints[maxLabel])
                                {
                                    outerBoundaryName = thisPatchName;
                                    outerBoundaryLabel = patchI;
                                }
                            }
                        }
                        Pout << "innerBoundaryName: " << innerBoundaryName << endl;
                        //Pout << "innerBoundaryLabel: " << innerBoundaryLabel << endl;
                        Pout << "outerBoundaryName: " << outerBoundaryName << endl;
                        //Pout << "outerBoundaryLabel: " << outerBoundaryLabel << endl;
                        // 8) From the remaining coil boundaries + windingDirectionAP find top and bottom boundary centers.
                        point topCenter(0,0,0);
                        point botCenter(0,0,0);
                        point topNormal(0,0,0);
                        point botNormal(0,0,0);
                        word topBoundaryName = "";
                        word botBoundaryName = "";
                        forAll(boundaryMesh,patchI)
                        {
                            word thisPatchName = boundaryMesh[patchI].name();
                            if (thisPatchName == terminalName || thisPatchName == groundTerminalName || thisPatchName == innerBoundaryName || thisPatchName == outerBoundaryName)
                                continue;
                            //Pout << "thisPatchName: " << thisPatchName << endl;
                            if (topBoundaryName == "")
                            {
                                topCenter = getCenter(regionName,thisPatchName);
                                topNormal = getNormal(regionName,thisPatchName);
                                topBoundaryName = thisPatchName;
                            }
                            else
                            {
                                botCenter = getCenter(regionName,thisPatchName);
                                botNormal = getNormal(regionName,thisPatchName);
                                botBoundaryName = thisPatchName;
                                break;
                            }
                        }
                        // 9) Define coil height as the distance between bot and top boundary centers.
                        scalar coilHeight = mag(topCenter-botCenter);
                        // 10 Define windingDirection as the vector from bot to top boundary center.
                        Foam::vector windingDirection = (topCenter-botCenter)/coilHeight;
                        //Info << "windingDirection: " << windingDirection << endl;
                        if ((topNormal & botNormal) < 0) botNormal*=-1;
                        Foam::vector coilNormal = 0.5*(topNormal + botNormal)/mag(0.5*(topNormal + botNormal));
                        if ((coilNormal & windingDirection) < 0) coilNormal*=-1;
                        Info << "coilNormal: " << coilNormal << endl;
                        if ( (windingDirection & windingDirectionAP) < 0)
                        {
                            //reverse
                            windingDirection *= -1;
                            point topCenterTemp = topCenter;
                            word topBoundaryNameTemp = topBoundaryName;
                            topCenter = botCenter;
                            topBoundaryName = botBoundaryName;
                            botCenter = topCenterTemp;
                            botBoundaryName = topBoundaryNameTemp;
                        }
                        windingDirection *= coilHeight;
                        Info << "windingDirection: " << windingDirection << endl;
                        Info << "topBoundaryName: " << topBoundaryName << endl;
                        Info << "topCenter: " << topCenter << endl;
                        Info << "botBoundaryName: " << botBoundaryName << endl;
                        Pout << "botCenter: " << botCenter << endl;
                        // 11) Construct line perpendicular to windingDirection and terminalDirection, which originates in present cell point.
                        // (Cross product between windingDirection and coilCenter to terminalCenter vector, terminalDirection,
                        // should be a safe line direction.)
                        Foam::vector search_vector = (coilNormal ^ terminalDirection)/mag(coilNormal ^ terminalDirection);//v1
                        //if (mag(search_vector) < SMALL) continue;
                        const word coilCoreName = terminalName+"_core";
                        Info << "search_vector: " << search_vector << endl;
                        // Inner boundary mesh
                        const fvPatch& inner_patch = boundaryMesh[innerBoundaryLabel];
                        // Outer boundary mesh
                        const fvPatch& outer_patch = boundaryMesh[outerBoundaryLabel];
                        forAll(names_, j)
                        {
                            const word& otherRegionName = names_[j].first();
                            if (otherRegionName == regionName)
                                continue;
                            if (!mesh(otherRegionName).objectRegistry::foundObject<volScalarField>(coilCoreName))
                            {
                                volScalarField* fPtr
                                (
                                    new volScalarField
                                    (
                                        IOobject
                                        (
                                            coilCoreName,
                                            mesh(otherRegionName).time().name(),
                                            mesh(otherRegionName),
                                            IOobject::NO_READ,
                                            IOobject::AUTO_WRITE
                                        ),
                                        mesh(otherRegionName),
                                        dimensionedScalar("",dimless,0)
                                    )
                                );

                                // Transfer ownership of this object to the objectRegistry
                                fPtr->store(fPtr);
                            }
                            volScalarField& coilCoreIds = mesh(otherRegionName).objectRegistry::lookupObjectRef<volScalarField>(coilCoreName);

                            vectorField regionCells = mesh(otherRegionName).C();
                            //Info << "Continuing for region: " << otherRegionName << endl;
                            int counter = 0;
                            forAll(regionCells,cellI)
                            {
                                const point test_point = regionCells[cellI];//p1
                                // 12) For each mesh surface on coil inner and outer boundary get center position, face_center,
                                // and using cell center coordinates as test_point and check if search_vector intersects a given face.
                                int inner_hit_count = surfaceHitCount_(test_point,search_vector,inner_patch);
                                int outer_hit_count = surfaceHitCount_(test_point,search_vector,outer_patch);
                                // 13) If line goes through exactly 1 inner and 1 outer boundary of coil, then cell is inside coil.
                                // TODO: Will count as a false negative, if intersection is exactly on a line segment or a vertex
                                // Could additionally use meshSearch.intersection() if inner_hit_count == 2 || outer_hit_count == 2
                                if (inner_hit_count == 1 && outer_hit_count == 1)
                                {
                                    coilCoreIds[cellI] = 1;
                                    counter++;
                                }
                            }
                            /*point testPoint(0,-0.15,-0.025);
                            int inner_hit_count = surfaceHitCount_(testPoint,search_vector,innerPatch);
                            int outer_hit_count = surfaceHitCount_(testPoint,search_vector,outerPatch);*/
                            Pout << "Region: " << otherRegionName << " cells " << regionCells.size() << " hits: " << counter << endl;
                            coilCoreIds.write();
                            // Write found properties to the dictionary of this region
                            //dictionary wireProperties;
                            //wireProperties.add("windingDirection", windingDirection, true);
                            //electroDict.add("wireProperties", wireProperties, true); //true flag merges dictionaries
                            electroDict.add("windingDirection", windingDirection, true); //true flag merges dictionaries
                            electroDict.regIOobject::write();
                            // Note: Electro dict values are rounded up to 6 significant digits after rewriting dict.
                        }
                    }

                    //Check if initialized
                    if (!electroDict.found("windingDirection"))//isDict("wireProperties"))
                    {
                        FatalIOError << "Failed to get windingDirection from " << regionName << " dictionary!\n"
                        << exit(FatalIOError);
                    }
                    forAll(names_, j)
                    {
                        const word& otherRegionName = names_[j].first();
                        if (otherRegionName == regionName)
                            continue;
                        const word coilCoreName = terminalName+"_core";
                        if (!mesh(otherRegionName).objectRegistry::foundObject<volScalarField>(coilCoreName))
                        {
                            volScalarField* fPtr
                            (
                                new volScalarField
                                (
                                    IOobject
                                    (
                                        coilCoreName,
                                        mesh(otherRegionName).time().name(),
                                        mesh(otherRegionName),
                                        IOobject::MUST_READ,
                                        IOobject::AUTO_WRITE
                                    ),
                                    mesh(otherRegionName)
                                )
                            );

                            // Transfer ownership of this object to the objectRegistry
                            fPtr->store(fPtr);
                        }
                    }
                    const scalar target_value =
                    dimensionedScalar(
                        "voltage",
                        dimMass*dimLength*dimLength/(dimTime*dimTime*dimTime*dimCurrent),
                        electroDict
                    ).value();
                    const scalar value_multiplier =
                    dimensionedScalar(
                        "voltageMultiplier",
                        dimless,
                        electrPtr->electro
                    ).value();// Number of turns / fill factor
                    const Foam::vector winding_direction = electroDict.lookup<Foam::vector>("windingDirection");
                    coilParameters_ terminalControls;
                    terminalControls.terminal_name = terminalName;
                    terminalControls.region_name = regionName;
                    terminalControls.target_value = target_value;
                    terminalControls.target_phase = target_phase;
                    terminalControls.value_multiplier = value_multiplier;
                    terminalControls.winding_direction = winding_direction;
                    voltageControls_[terminalName] = terminalControls;
                    voltageNeedsUpdate_[terminalName] = true;
                    scalar oldVoltageValue = readControlValue_("coilVoltages/"+terminalName);
                    scalar oldVoltagePhase = readControlValue_("coilPhases/"+terminalName);
                    inducedVoltageValue_[terminalName] = oldVoltageValue - target_value;
                    inducedVoltagePhase_[terminalName] = oldVoltagePhase - target_phase;
                    Info << "inducedVoltageValue_: " << inducedVoltageValue_[terminalName] << endl;
                    Info << "inducedVoltagePhase_: " << inducedVoltagePhase_[terminalName] << endl;
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
    if (isFluid(regionName) || isSolid(regionName) || isElectric(regionName) || isMagnetic(regionName))
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
Foam::scalar Foam::conductingRegionSolvers::getPotErefValue_(const word terminalName, bool imaginary)
{
    scalar refValue = 0;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            electroBase* electrPtr = getElectroBasePtr_(regionName);
            if (electrPtr->hasBoundary(terminalName))
            {
                refValue = electrPtr->getPotErefValue(terminalName, imaginary);
            }
        }
    }
    return refValue;
}
void Foam::conductingRegionSolvers::updatePotErefValue_(const word terminalName, scalar newValue, bool imaginary)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            getElectroBasePtr_(regionName)->updatePotErefValue(terminalName, newValue, imaginary);
        }
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

void Foam::conductingRegionSolvers::checkIfAnyVoltageControl_()
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getElectroBasePtr_(regionName))
        {
            electroBase* electrPtr = getElectroBasePtr_(regionName);
            const word feedbackType = electrPtr->electro.lookupOrDefault<word>("feedbackControl","");
            if (feedbackType == "voltage")
            {
                hasVoltageControl_ = true;
                break;
            }
        }
    }
}

int Foam::conductingRegionSolvers::surfaceHitCount_(const point test_point,const Foam::vector search_vector,const fvPatch& patch)
{
    const vectorField centers = patch.Cf();
    const faceList& faces = patch.patch().localFaces();
    const pointField& points = patch.patch().localPoints();
    int outer_hit_count = 0;

    forAll(faces, faceI)
    {
        const face& test_face = faces[faceI];
        const point face_center = centers[faceI];//p2

        // Defensive checks
        if (test_face.size() < 3) continue;//Not a face
        if (max(test_face) >= points.size()) continue;//Out of bounds check
        if (((test_point-face_center)&search_vector) < 0) continue;//wrong direction

        // Triangulate face and check for a hit
        point p0 = points[test_face[0]];
        for (label i = 1; i < test_face.size() - 1; ++i)
        {
            point p1 = points[test_face[i]];
            point p2 = points[test_face[i+1]];
        
            Foam::triPointRef tri(p0, p1, p2);
        
            pointHit hit = tri.ray(test_point, search_vector);
        
            if (hit.hit())
            {
                outer_hit_count++;
                break;
            }
        }
    }
    return outer_hit_count;
}

Foam::volVectorField Foam::conductingRegionSolvers::getJdirection_(word regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName) && getElectroBasePtr_(regionName)->currentReferenceSet(imaginary))
    {
        volVectorField regionField = getElectroBasePtr_(regionName)->getJref(imaginary);
        volScalarField magnitudeField = mag(regionField);
        if (max(magnitudeField).value() == 0)
        {
            FatalIOError << "Region" << regionName 
            << "reference field for J" << (isElectroHarmonic() ? (imaginary ? "im" : "re") : "") << " is equal to zero!\n"
            << exit(FatalIOError);
        }
        volVectorField directionField = regionField/magnitudeField;
        return directionField;
    }
    else
    {
        FatalIOError << "Failed to get J from region " << regionName << "!\n"
        << exit(FatalIOError);
    }
}

Foam::scalar Foam::conductingRegionSolvers::getInductionSum_(word regionName, const coilParameters_ terminalControls, bool imaginary)
{
    Foam::vector winding_direction = terminalControls.winding_direction;
    const scalarField volume = mesh(regionName).V();
    bool coil_region = terminalControls.region_name == regionName;
    const word coilCoreName = terminalControls.terminal_name+"_core";
    const scalarField coreFilter = coil_region ?
        scalarField(volume.size(), 1.0) :
        scalarField(mesh(regionName).objectRegistry::lookupObject<volScalarField>(coilCoreName).internalField());
    scalar sumI = 0;
    if (isElectroHarmonic())
    {
        const scalarField Bnow = getElectro(regionName).B(!imaginary).internalField() & winding_direction;
        //gSum() (also gAverage()) is multi-threading-safe way to sum over all grid points.
        sumI = gSum(volume*Bnow*coreFilter);//Find volume integral
    }
    else
    {
        const scalarField Bnow = getElectro(regionName).B(imaginary).internalField() & winding_direction;
        const scalarField Bold = regionOldB_[regionName] & winding_direction;
        //gSum() (also gAverage()) is multi-threading-safe way to sum over all grid points.
        sumI = gSum(volume*Bnow*coreFilter)-gSum(volume*Bold*coreFilter);//Find volume integral
    }
    //Info << "region: " << regionName << " volume: " << gSum(volume) << endl;
    //Info << "sumI: " << sumI << endl;
    return sumI;
}

Foam::scalar Foam::conductingRegionSolvers::getCurrentSum_(word regionName,bool imaginary)//volVectorField JGlobal,
{
    const volVectorField Jdirection = getJdirection_(regionName,imaginary);//Get reference
    const volVectorField JRegion = getElectro(regionName).J(imaginary);//Use copy of internal J for rewriting (has respective regions' mesh)
    const scalarField Jprojection(JRegion & Jdirection);//Project to a reference
    const scalarField volume = mesh(regionName).V();
    //gSum() (also gAverage()) is multi-threading-safe way to sum over all grid points.
    scalar sumJ = gSum(volume*Jprojection)/max(gSum(volume), vSmall);//Find volume average
    return sumJ;
}

Foam::scalar Foam::conductingRegionSolvers::getTransientPhaseShift_(scalar argument, scalar step_shift)
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
    scalar time_position = runTime_.userTimeValue() - step_shift*tStepElectro_;
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

bool Foam::conductingRegionSolvers::fileExists_(const word fileName)
{
    bool file_exists = false;
    std::ifstream inputFile(fileName, std::ios::in);
    if (inputFile.is_open())
    {
        file_exists = true;
        inputFile.close();
    }
    return file_exists;
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

Foam::scalar Foam::conductingRegionSolvers::getArea(const word regionName, const word patchName)
{
    const fvBoundaryMesh& boundaryMesh = mesh(regionName).boundary();
    forAll(boundaryMesh, patchi)
    {
        const fvPatch& thisPatch = boundaryMesh[patchi];
        if (thisPatch.name() == patchName)
        {
            const scalarField surface = thisPatch.magSf();
            scalar totVol = gSum(surface);
            return totVol;
        }
    }
    FatalIOError << "Region " << regionName << " doesn't have a boundary named " << patchName << "!\n"
    << exit(FatalIOError);
    return -1;
}

Foam::vector Foam::conductingRegionSolvers::getCenter(const word regionName)
{
    const scalarField volume = mesh(regionName).V();
    const vectorField cellPoints = mesh(regionName).C();
    const scalarField points_x = cellPoints & Foam::vector(1,0,0);
    const scalarField points_y = cellPoints & Foam::vector(0,1,0);
    const scalarField points_z = cellPoints & Foam::vector(0,0,1);
    scalar totVol = gSum(volume);
    scalar totVolx = gSum(volume*points_x)/totVol;
    scalar totVoly = gSum(volume*points_y)/totVol;
    scalar totVolz = gSum(volume*points_z)/totVol;
    Foam::vector coilCenter(totVolx,totVoly,totVolz);
    return coilCenter;
}

Foam::vector Foam::conductingRegionSolvers::getCenter(const word regionName, const word patchName)
{
    //Pout << "In get center" << endl;
    const fvBoundaryMesh& boundaryMesh = mesh(regionName).boundary();
    forAll(boundaryMesh, patchi)
    {
        const fvPatch& thisPatch = boundaryMesh[patchi];
        if (thisPatch.name() == patchName)
        {
            const scalarField surface = thisPatch.magSf();
            const vectorField points = thisPatch.Cf();
            const scalarField points_x = points & Foam::vector(1,0,0);
            const scalarField points_y = points & Foam::vector(0,1,0);
            const scalarField points_z = points & Foam::vector(0,0,1);
            scalar totVol = gSum(surface);
            scalar totVolx = gSum(surface*points_x)/totVol;
            scalar totVoly = gSum(surface*points_y)/totVol;
            scalar totVolz = gSum(surface*points_z)/totVol;
            Foam::vector centerPosition(totVolx,totVoly,totVolz);
            return centerPosition;
        }
    }
    FatalIOError << "Region " << regionName << " doesn't have a boundary named " << patchName << "!\n"
    << exit(FatalIOError);
    return Foam::vector(0,0,0);
}

Foam::vector Foam::conductingRegionSolvers::getNormal(const word regionName, const word patchName)
{
    const fvBoundaryMesh& boundaryMesh = mesh(regionName).boundary();
    forAll(boundaryMesh, patchi)
    {
        const fvPatch& thisPatch = boundaryMesh[patchi];
        if (thisPatch.name() == patchName)
        {
            Pout << "thisPatch.size() " << thisPatch.size() << endl;
            const scalarField surface = thisPatch.magSf();
            const vectorField points = thisPatch.nf();
            const scalarField points_x = points & Foam::vector(1,0,0);
            const scalarField points_y = points & Foam::vector(0,1,0);
            const scalarField points_z = points & Foam::vector(0,0,1);
            scalar totVol = gSum(surface);
            scalar totVolx = gSum(surface*points_x)/totVol;
            scalar totVoly = gSum(surface*points_y)/totVol;
            scalar totVolz = gSum(surface*points_z)/totVol;
            Foam::vector averageNormal(totVolx,totVoly,totVolz);
            return averageNormal;
        }
    }
    FatalIOError << "Region " << regionName << " doesn't have a boundary named " << patchName << "!\n"
    << exit(FatalIOError);
    return Foam::vector(0,0,0);
}

bool Foam::conductingRegionSolvers::checkIfAnyVelocityRegions()
{
    bool hasVelocityRegions = false;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (isFluid(regionName))
        {
            hasVelocityRegions = true;
            break;
        }
    }
    return hasVelocityRegions;
}

void Foam::conductingRegionSolvers::updateFeedbackControl()//volVectorField& JreGlobal,volVectorField& JimGlobal
{
    // Elmer should be called only at writeTime with adjustableTime option enabled.
    // Check if for some reason that is not so.
    //if ( adjustableRunTime_ && !runTime_.writeTime() ) {return;}
    // TODO: don't increase time-step size if needsUpdate
    // OR run extra initial iterations to get to correct values.
    if (!isElectroHarmonic())
    {
        integrationCounter_++;
    }
    // Calculate total induced voltage from Faraday's law of induction
    // value_multiplier = number_of_turns/fill_factor
    // induced_voltage = -value_multiplier*Integral((dB/dt)*d(area_perpendicular))
    // induced_voltage is applied for the whole coil so assuming == const.
    // Additionally, calcuate integral average over coil height to get the effective voltage for coil terminal.
    // coil_height = mag(winding_direction)
    // induced_voltage *coil_height = -value_multiplier*Integral((dB/dt)&(winding_direction/coil_height)*d(volume))
    // Coil has finite width, so let's calculate average over volume inside coil (inner_core)
    // and total volume (inner_core + coil_volume).
    // Let's say, integralCore = Integral((dB/dt)*d(core_volume))
    // and integralCoil = Integral((dB/dt)*d(coil_volume))
    // We get:
    // induced_voltage = -(integralCore+0.5*integralCoil)*value_multiplier/coil_height
    // Update voltage on the coil terminal as the initially applied voltage plus the induced voltage.
    // Continue until change of induced voltage between two consecutive steps is sufficiently small.
    for (auto element : terminalToRegions_)
    {
        // Coil regions
        const word terminalName = element.first;
        Info << "terminalName: " << terminalName << endl;
        if (voltageControls_.find(terminalName) != voltageControls_.end() && voltageNeedsUpdate_[terminalName])
        {
            // In harmonic case
            // real part of induction is stored in inducedVoltageValue_
            // and imaginary part in inducedVoltagePhase_
            // In transient case
            // Integration sum of deltaB * winding_direction is stored in inducedVoltageValue_
            // The induction component of the present time step is stored in inducedVoltagePhase_

            const coilParameters_ terminalControls = voltageControls_[terminalName];
            for (word regionName : element.second)
            {
                // Reseting also for transient case, because used for storing momentary induction.
                inducedVoltagePhase_[terminalName] = 0;
                if (isElectroHarmonic())
                {
                    inducedVoltageValue_[terminalName] = 0;
                }
                forAll(names_, i)
                {
                    const word& otherRegionName = names_[i].first();
                    // Coefficients according to
                    // induced_voltage = -(integralCore+0.5*integralCoil)*value_multiplier/coil_height
                    scalar coeff = -terminalControls.value_multiplier/mag(terminalControls.winding_direction);
                    if (otherRegionName == regionName)//coil/wire region
                    {
                        coeff *= 0.5;
                    }

                    if (isElectroHarmonic())
                    {
                        // Harmonic case
                        //dBdt = omega*i*B
                        //dBdt Re = -omega*Bim
                        //dBdt Im = omega*Bre

                        // Real part
                        scalar time_derivative = -2*PI*frequency_;
                        const scalar magneticIntegralRe = getInductionSum_(otherRegionName,terminalControls);
                        inducedVoltageValue_[terminalName] += coeff*time_derivative*magneticIntegralRe;

                        // Imaginary part
                        time_derivative = 2*PI*frequency_;
                        const scalar magneticIntegralIm = getInductionSum_(otherRegionName,terminalControls,true);
                        inducedVoltagePhase_[terminalName] += coeff*time_derivative*magneticIntegralIm;
                    }
                    else
                    {
                        // Transient case
                        // Integrate over time to get magnitude
                        // d/dt for DB/dt is
                        scalar time_derivative = 1.0/tStepElectro_;
                        const scalar induction_component =
                        coeff*time_derivative*getInductionSum_(otherRegionName,terminalControls);
                        //inducedVoltageValue_[terminalName] += abs(induction_component)*tStepElectro_;//Integration component
                        inducedVoltagePhase_[terminalName] += induction_component;//Momentary induction component
                        const vectorField oldB = getElectro(otherRegionName).B().internalField();
                        regionOldB_[otherRegionName]=oldB;
                    }
                }
                if (!isElectroHarmonic())
                {
                    inducedVoltageValue_[terminalName] += abs(inducedVoltagePhase_[terminalName])*tStepElectro_;//Integration component
                }
            }
            //Convert to value and phase
            if (isElectroHarmonic() || integrationCounter_ % integrationSteps_ == 0)
            {
                // Harmonic case or half period of transient case
                Info << "Starting induced voltage update " << endl;
                Info << "Updating voltage for terminal " << terminalName << endl;
                const scalar oldVoltageValue = readControlValue_("coilVoltages/"+terminalName);
                const scalar oldVoltagePhase = readControlValue_("coilPhases/"+terminalName);
                const scalar setVoltageRe = terminalControls.target_value*Foam::cos(terminalControls.target_phase*PI/180.0);
                const scalar setVoltageIm = terminalControls.target_value*Foam::sin(terminalControls.target_phase*PI/180.0);
                const scalar oldVoltageRe = oldVoltageValue*Foam::cos(oldVoltagePhase*PI/180.0);
                const scalar oldVoltageIm = oldVoltageValue*Foam::sin(oldVoltagePhase*PI/180.0);
                // Calculate previous induction
                const scalar oldInducedVoltageRe = oldVoltageRe - setVoltageRe;
                const scalar oldInducedVoltageIm = oldVoltageIm - setVoltageIm;

                // Induced values
                scalar newInducedVoltageRe = 0;
                scalar newInducedVoltageIm = 0;
                if (isElectroHarmonic())
                {
                    newInducedVoltageRe = inducedVoltageValue_[terminalName];
                    newInducedVoltageIm = inducedVoltagePhase_[terminalName];
                }
                else
                {
                    // Calculate induced voltage magnitude and phase from integral and latest time value
                    // Integral (A*sin(2*pi*f*t)) from 0 to 1/(2*f) = A/(pi*f)
                    inducedVoltageValue_[terminalName] *= PI*frequency_;
                    scalar inducedVoltageValue = inducedVoltageValue_[terminalName];
                    scalar inducedVoltagePhase = getTransientPhaseShift_(
                        inducedVoltagePhase_[terminalName]//Induction value at latest time step
                        /inducedVoltageValue,//Integral value over integration time
                    0.5)*180/PI;//TODO: Check if using half time-step shift works reliably
                    
                    newInducedVoltageRe = inducedVoltageValue*Foam::cos(inducedVoltagePhase*PI/180.0);
                    newInducedVoltageIm = inducedVoltageValue*Foam::sin(inducedVoltagePhase*PI/180.0);
                    /*Info << "integral induction: " << inducedVoltageValue << endl;
                    Info << "momentary induction: " << inducedVoltagePhase_[terminalName] << endl;
                    Info << "division : " << inducedVoltagePhase_[terminalName]/inducedVoltageValue << endl;
                    Info << "Induced "  << endl;
                    Info << "inducedVoltageValue: " << inducedVoltageValue << endl;
                    Info << "inducedVoltagePhase: " << inducedVoltagePhase << endl;*/
                    // Reset for next integration
                    inducedVoltageValue_[terminalName] = 0;
                    inducedVoltagePhase_[terminalName] = 0;
                }
                //Info << "inducedVoltageRe: " << newInducedVoltageRe << endl;
                //Info << "inducedVoltageIm: " << newInducedVoltageIm << endl;

                // Calculate induction difference
                scalar deltaInducedRe = newInducedVoltageRe-oldInducedVoltageRe;
                scalar deltaInducedIm = newInducedVoltageIm-oldInducedVoltageIm;
                //Info << "deltaInducedRe: " << deltaInducedRe << endl;
                //Info << "deltaInducedIm: " << deltaInducedIm << endl;

                // Apply relaxation to induced value difference to improve control stability.
                deltaInducedRe *= induced_voltage_relaxation_factor_;
                deltaInducedIm *= induced_voltage_relaxation_factor_;
                scalar newVoltageRe = setVoltageRe// Adding preset voltage value
                + oldInducedVoltageRe// Adding previous induced voltage value
                + deltaInducedRe// Adding voltage difference value induced by previous imaginary value
                //Account for backreaction on next step
                + deltaInducedIm//New imaginary induced by previous real
                *(
                    (abs(oldVoltageIm) > SMALL)
                    ?
                    (newInducedVoltageRe/oldVoltageIm)//Weight of how much imaginary part induces real
                    :
                    //Next best guess if oldVoltageIm unavailable
                    (sqrt(pow(newInducedVoltageRe,2)+pow(newInducedVoltageIm,2))/oldVoltageValue)
                )*induced_voltage_relaxation_factor_;//This is what will be induced on the next step, so relax that as well.
                scalar newVoltageIm = setVoltageIm// Adding preset voltage value
                + oldInducedVoltageIm// Adding previous induced voltage value
                + deltaInducedIm// Adding voltage difference value induced by previous real value
                //Account for backreaction on next step
                - deltaInducedRe//New real induced by previous imaginary
                *(
                    abs(oldVoltageRe) > SMALL
                    ?
                    (newInducedVoltageIm/oldVoltageRe)//Weight of how much real part induces imaginary
                    :
                    //Next best guess if oldVoltageRe unavailable
                    (sqrt(pow(newInducedVoltageRe,2)+pow(newInducedVoltageIm,2))/oldVoltageValue)
                )*induced_voltage_relaxation_factor_;//This is what will be induced on the next step, so relax that as well.

                scalar newVoltageValue = sqrt(pow(newVoltageRe,2)+pow(newVoltageIm,2));
                scalar newVoltagePhase = atan2(newVoltageIm,newVoltageRe)*180/PI;
                //Info << "newVoltageRe: " << newVoltageRe << endl;
                //Info << "newVoltageIm: " << newVoltageIm << endl;

                scalar errorValue = abs(newVoltageValue - oldVoltageValue)/max(abs(newVoltageValue),vSmall);
                scalar errorPhase = abs(newVoltagePhase - oldVoltagePhase)/180.0;

                Info << "Previous voltage value: " << oldVoltageValue << "; phase: " << oldVoltagePhase << endl
                << "New voltage value: " << newVoltageValue << "; phase: " << newVoltagePhase << endl
                << "Value error: " << errorValue << "; phase error: " << errorPhase << endl;

                const bool quit_requirement =
                (isElectroHarmonic()
                ?
                initial_iter_ > 1//Has one extra initial iteration for harmonic case, because first error will be zero.
                :
                initial_iter_ > 1 && wait_iter_ > waitInterval);//First go through initialization.

                if (errorValue < target_value_error_ && errorPhase < target_phase_error_ && quit_requirement)
                {
                    voltageNeedsUpdate_[terminalName] = false;
                    Info << "Voltage requirements met!" << endl;
                    continue;
                }
                const int max_initial_iter = isElectroHarmonic() ? 1 : waitInterval+1;
                if (initial_iter_ <= max_initial_iter)
                {
                    initial_iter_++;
                }
                // Write new voltage values to file
                writeControlValue_("coilVoltages/"+terminalName,newVoltageValue);
                writeControlValue_("coilPhases/"+terminalName,newVoltagePhase);
            }
            
            if (!isElectroHarmonic())
            {
                /*if (integrationCounter_ % integrationSteps_ != 0)
                {
                    // Intermediate step between half-periods
                    Info << "integral: " << inducedVoltageValue_[terminalName] << endl;
                    Info << "momentary induction: " << inducedVoltagePhase_[terminalName] << endl;
                }*/
                if (wait_iter_ <= waitInterval)
                {
                    wait_iter_++;
                }
            }
        }
    }
    // Calculate the integral current density and compare it to target value
    // Update voltage on the wire terminal through the control of feedback loop
    // Continue until error is sufficiently small.
    for (auto element : terminalToRegions_)
    {
        const word terminalName = element.first;
        if (feedbackControllers_.find(terminalName) != feedbackControllers_.end() && feedbackControllers_[terminalName].getControlType() == "current")
        {
            scalar avgJre = 0;
            scalar avgJim = 0;
            for (word regionName : element.second)
            {
                avgJre += getCurrentSum_(regionName);
                if (isElectroHarmonic())
                {
                    avgJim += getCurrentSum_(regionName,true);
                }
            }
            
            if (!isElectroHarmonic())
            {
                if (integrationCounter_ % integrationSteps_ != 0)
                {
                    integralCurrent_[terminalName] += abs(avgJre)*tStepElectro_;
                    return;
                }
                integralCurrent_[terminalName] += abs(avgJre)*tStepElectro_;
                // Integral (A*sin(2*pi*f*t)) from 0 to 1/(2*f) = A/(pi*f)
                integralCurrent_[terminalName] *= PI*frequency_;
            }
            
            const scalar present_value = isElectroHarmonic() ? sqrt(pow(avgJre,2)+pow(avgJim,2))
            : integralCurrent_[terminalName];
            const scalar present_phase = isElectroHarmonic() ?
            atan2(avgJim,avgJre)*180/PI :
            getTransientPhaseShift_(avgJre/integralCurrent_[terminalName],0.5)*180/PI;//This is at midpoint of avgJre
            integralCurrent_[terminalName] = 0;
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
                writeControlValue_("coilVoltages/"+terminalName,control_values.first());
                writeControlValue_("coilPhases/"+terminalName,control_values.second());
                
                // Current value raises very slowly for high frequencies, because the initial estimate is wrong by several orders of magnitude.
                // Updating coefficient later helps with that.
                Pair<scalar> old_proportional_coeff = feedbackControllers_[terminalName].getProportionalCoefficient();
                scalar target_value = feedbackControllers_[terminalName].getReference().first();
                scalar old_value_coeff = old_proportional_coeff.first();
                scalar new_value_coeff = -0.5*control_values.first()/target_value;
                bool needs_coeff_update = abs(new_value_coeff) > 10*old_value_coeff;
                if (needs_coeff_update)
                {
                    Pair<scalar> new_proportional_coeff(new_value_coeff,old_proportional_coeff.second());
                    feedbackControllers_[terminalName].updateProportionalCoefficient(new_proportional_coeff);
                }
            }
            else
            {
                Info << "Calculated current for terminal " << terminalName  << " is within acceptable errors." << endl;
            }
        }
    }
}

bool Foam::conductingRegionSolvers::controllersNeedUpdate()
{
    for (auto element : terminalToRegions_)
    {
        const word terminalName = element.first;
        // Current controllers
        if (feedbackControllers_.find(terminalName) != feedbackControllers_.end() && feedbackControllers_[terminalName].needsUpdate())
        {return true;}
        // Voltage controllers
        if (voltageControls_.find(terminalName) != voltageControls_.end() && voltageNeedsUpdate_[terminalName])
        {return true;}
    }
    return false;
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
/*bool Foam::conductingRegionSolvers::hasElectricSources()
{
    return hasElectricSources_;
}*/
//Get J from electromagnetic source region
/*Foam::volVectorField Foam::conductingRegionSolvers::getSourceJ(word regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName) && getElectroBasePtr_(regionName)->electro.getRegionRole() == "coil")
    {
        electroBase* electrPtr = getElectroBasePtr_(regionName);
        volVectorField regionField = isSource(regionName) ? electrPtr->getJ(imaginary) : electrPtr->getJref(imaginary);
        if (!isElectroHarmonic())
        {
                regionField*=conductorValues_[regionName]*Foam::sin(2*PI*frequency_*runTime_.userTimeValue()
                //TODO: Check if sign is consistent everywhere
                -conductorPhases_[regionName]*PI/180.0);
        }
        return regionField;
    }
    else
    {
        FatalIOError << "Failed to get J from region " << regionName << "!\n"
        << exit(FatalIOError);
    }
}*/
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
//Assigns fluid and solid region values from global field to a specified region field
void Foam::conductingRegionSolvers::setJRefToRegion(volVectorField& globalField, const word& regionName, bool imaginary)
{
    if (getElectroBasePtr_(regionName))
    {
        electroBase* electrPtr = getElectroBasePtr_(regionName);
        if (electrPtr->currentReferenceSet(imaginary)) return;// Set only once
        volVectorField& regionField = electrPtr->getJref(imaginary);
        vectorGlobalToField_(globalField,regionField,regionName);
        // Update boundary values based on boundary conditions.
        regionField.correctBoundaryConditions();
        electrPtr->markCurrentReferenceAsSet(imaginary);
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

bool Foam::conductingRegionSolvers::isMagnetic(const word regionName)
{
    return names_[regionIdx_[regionName]].second() == magneticSolverName_;
}

bool Foam::conductingRegionSolvers::needsReference(const word regionName)
{
    const electromagneticModel& electro = getElectro(regionName);
    const word regionRole = electro.getRegionRole();
    const word feedbackType = electro.lookupOrDefault<word>("feedbackControl","");
    const bool regionNeedsReference = isElectric(regionName) && (regionRole == "wire") && (feedbackType == "current");
    if (!regionNeedsReference) return false;
    const bool referenceNotSet = isElectroHarmonic()
    ?
    !getElectroBasePtr_(regionName)->currentReferenceSet() || !getElectroBasePtr_(regionName)->currentReferenceSet(true)
    :
    !getElectroBasePtr_(regionName)->currentReferenceSet();
    return referenceNotSet;
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
                scalar maxMagneticReynoldsDifference = mu_0 * characteristicSizes_[i] *
                    gMax((getElectro(regionName).sigmaInv()*mag(U_old-U))());
                scalar maxRelativeVelocityDifference = gMax((mag(U_old-U)/(gAverage((mag(U))())+smallU))());
                maxRemDiff_local = max(
                    maxMagneticReynoldsDifference,
                    maxRemDiff_local);        

                maxRelDiff_local = max(
                    maxRelativeVelocityDifference,
                    maxRemDiff_local);
            }
        }

        if((maxRelDiff_local>maxRelDiff_ || maxRelDiff_<SMALL) && maxRelDiff_+SMALL<=1.0) {
            doUpdate = true;
        }
        else if(maxRemDiff_local>maxRemDiff_ && maxRelDiff_-SMALL<=1.0) {
            doUpdate = true;
        }
        // Include update check for controllers
        for (auto element : terminalToRegions_)
        {
            const word terminalName = element.first;
            Info << "terminalName: " << terminalName << endl;
            if (voltageControls_.find(terminalName) != voltageControls_.end() && voltageNeedsUpdate_[terminalName])
            {doUpdate = true;}
            if (feedbackControllers_.find(terminalName) != feedbackControllers_.end()  && feedbackControllers_[terminalName].needsUpdate())
            {doUpdate = true;}
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
