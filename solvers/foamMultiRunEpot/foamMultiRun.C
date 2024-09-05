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

Application
    foamMultiRunEpot is a modified foamMultiRun solver, where the 
    modification is based on EOF-Library solver mhdVxBPimpleFoam. Additional 
    modification was made to update electrical currents locally in OpenFOAM,
    for currents, which doesn't disturb external magnetic field. More precisely,
    while the change in magnetic Reynolds number doesn't exceed the provided value.

Description
    Solver for steady or transient electromagnetically forced fluid flow and 
    solid heat conduction, with conjugate heat transfer between regions, 
    buoyancy effects, turbulence and radiation modelling.
    Implements coupling with ElmerFEM solver.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "conductingRegionSolvers.H"
#include "pimpleMultiRegionControl.H"
#include "setDeltaT.H"

using namespace Foam;
#include "Elmer.H"
#include <fstream>
#include "fieldMapper.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // Create the region meshes and solvers
    conductingRegionSolvers solvers(runTime);

    // Create the outer PIMPLE loop and control structure
    pimpleMultiRegionControl pimple(runTime, solvers);
    
    // Set the initial time-step
    setDeltaT(runTime, solvers);

    // Read global mesh of all regions
    const fvMesh& meshGlobal = solvers.globalMesh();
    regionToGlobalCellId = solvers.regionToGlobalCellId;
    forAll(solvers.getNames(), i)
    {
        regionNames.append(solvers.getNames()[i].first());
    }

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //
    // EOF-Library supports using single mesh of communication either for sending or receiving
    //
    Info<< "*** Preparing Elmer communications for sending" << nl << endl;
    Elmer<fvMesh> sending(
            meshGlobal, //mesh
            1, // 1=send, -1=receive
            0, // 1=initialize, 0=w/o init
            1 // 1=multiregion, 0=exports O2E files
        );

    Info<< "*** Preparing Elmer communications for receiving" << nl << endl;
    Elmer<fvMesh> receiving(
            meshGlobal, //mesh
            -1, // 1=send, -1=receive
            0, // 1=initialize, 0=w/o init
            1 // 1=multiregion, 0=exports O2E files
        );

    double elmerClock = runTime.clockTimeIncrement();

    Info<< nl << "Initializing electromagnetic solver\n" << endl;

    // Create file for logging simulation times whenever Elmer is called
    string elmerTimesFileName = "postProcessing/elmerTimes.log";
    int elmer_status = 1; // 1=ok, 0=lastIter, -1=error
    bool initialize_elmer = true;
    #include "runElmerUpdate.H"
    initialize_elmer = false;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Write initial values
    #include "writeIntegrals.H"
    bool lastTimeStep = false;

    double OFClock = 0;

    Info<< nl << "Starting time loop\n" << endl;
    while (pimple.run(runTime) || lastTimeStep)
    {
        forAll(solvers, i)
        {
            solvers[i].preSolve();
        }

        solvers.setGlobalPrefix();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solvers);
        // Adjust time step so that last step is at end time.
        if (runTime.userTimeValue() + runTime.deltaTValue() > runTime.endTime().value())
        {
            const scalar lastDeltaT = runTime.endTime().value() - runTime.userTimeValue();
            runTime.setDeltaT(lastDeltaT);
            Info<< "Adjusting time step to match end time." << nl << endl;
            Info<< "deltaT = " << runTime.deltaTValue()  << nl << endl;
            lastTimeStep = true;
        }
        // Update paths for cleanup
        forAll(regionNames, i)
        {
            regionPaths[i] = getFieldPaths(solvers.mesh(regionNames[i]));
        }

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Multi-region PIMPLE corrector loop
        while (pimple.loop())
        {
            forAll(solvers, i)
            {
                solvers[i].moveMesh();
            }

            forAll(solvers, i)
            {
                solvers[i].fvModels().correct();
            }

            forAll(solvers, i)
            {
                solvers[i].prePredictor();
            }

            forAll(solvers, i)
            {
                solvers[i].momentumPredictor();
            }

            while (pimple.correctEnergy())
            {
                forAll(solvers, i)
                {
                    solvers[i].thermophysicalPredictor();
                }
            }

            forAll(solvers, i)
            {
                solvers[i].pressureCorrector();
            }

            forAll(solvers, i)
            {
                solvers[i].postCorrector();
            }
        }

        forAll(solvers, i)
        {
            solvers[i].postSolve();
        }
        //Update liquid-solid phase fraction
        forAll(regionNames, i)
        {
            if (solidificationEnabled[i])
            {
                alpha1Region[i] = solvers.mesh(regionNames[i]).lookupObject<volScalarField>(solverSolidificationName);
            }
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Check whether we need to update electromagnetic stuff with Elmer
        // Update electromagnetics with Elmer if magnetic field is
        // significantly disturbed and needs to be updated.
        // Otherwise, update electromagnetics locally by calculating electric potential.
        bool doElmer = solvers.updateMagneticField();

        // Calculate electric potential if current density will not be updated
        if (!doElmer)
        {
            solvers.setCorrectElectromagnetics();
        }

        forAll(regionNames, i)
        {
            solvers.solveElectromagnetics(regionNames[i]);
            const volVectorField& JxBRegion = solvers.getElectro(regionNames[i]).JxB;
            solvers.vectorFieldToGlobal(
                JxBGlobal,
                JxBRegion,
                regionNames[i]
            );
            const volScalarField& JJsigmaRegion = solvers.getElectro(regionNames[i]).JJsigma;
            solvers.scalarFieldToGlobal(
                JJsigmaGlobal,
                JJsigmaRegion,
                regionNames[i]
            );
        }

        solvers.setGlobalPrefix();
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Write last time step even if not write time
        if(runTime.writeTime() || lastTimeStep)
        {
            runTime.writeNow();
            // Cleanup
            forAll(regionNames, i)
            {
                forAll(regionPaths[i], j)
                {
                    if (!keepField[Pair<word>(regionPaths[i][j].first(),regionNames[i])])
                    {
                        //Pout << "Deleting file " << regionPaths[i][j].first() <<" was " <<
                        fileHandler().rm(regionPaths[i][j].second());// << endl;
                    }
                }
            }
        }
        //write integral values for all time steps
        #include "writeIntegrals.H"
        OFClock = runTime.clockTimeIncrement();

        Info<< "ExecutionTime : " << "Hydrodynamics step = " << OFClock << " s"
            << " ; Electrodynamics step = " << elmerClock << " s"
            << " ; ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
			
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update electromagnetic stuff with Elmer

        if(doElmer && runTime.run())
        {
            #include "runElmerUpdate.H"
        }
        elmerClock = runTime.clockTimeIncrement();
        // If run loop exited just before end time, schedule one more iteration.
        if (
            // Loop has been stopped.
            !runTime.run() &&
            // End time has not been reached.
            runTime.userTimeValue() < runTime.endTime().value() &&
            // Next step reaches end time.
            runTime.userTimeValue() + runTime.deltaTValue() >= runTime.endTime().value()
            )
        {
            lastTimeStep = true;
        }
        // Exit after extra iteration.
        else if (lastTimeStep)
        {
            lastTimeStep = false;
        }
    }

    Info << "Final iteration: "
        << " ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

    //Final iter for Elmer
    elmer_status = 0; // 1=ok, 0=lastIter, -1=error
    #include "runElmerUpdate.H"

    int clockDays = std::floor(runTime.elapsedClockTime()/3600.0/24.0);
    int clockHours = std::floor(runTime.elapsedClockTime()/3600.0-clockDays*24.0);
    int clockMinutes = std::floor(runTime.elapsedClockTime()/60.0-clockHours*60.0-clockDays*60.0*24.0);
    int clockSeconds = std::floor(runTime.elapsedClockTime()-clockMinutes*60.0-clockHours*3600.0-clockDays*3600.0*24.0);

	Info<< "Calculation completed in "
		<< clockDays << " days "
		<< clockHours << " h "
		<< clockMinutes << " min "  << clockSeconds << " s" 
		<< nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
