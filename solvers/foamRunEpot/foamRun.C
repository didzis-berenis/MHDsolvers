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
    foamRunEpot is a modified foamRun solver, where the modification 
    is based on EOF-Library solver mhdVxBPimpleFoam. Additional modification 
	was made to update electrical currents in OpenFOAM, while the change in 
	magnetic Reynolds number doesn't exceed the provided value. This 
	modification was based on the epotFoam solver, which can be found
	in https://doi.org/10.13140/RG.2.2.12839.55201 (Chapter 4).

Description
    Solver for steady or transient buoyant, turbulent flow of compressible
    or incompressible fluids for electromagnetically forced and heated flows.

    Compile option ELMER_TIME == HARMONIC_TIME builds foamRunEpot solver,
    which assumes coupling with harmonic (time-averaged) ElmerFEM solver.

    Compile option ELMER_TIME == TRANSIENT_TIME builds foamRunEpotTransient
    solver, which assumes coupling with transient ElmerFEM solver.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "conductingRegionSolver.H"
#include "pimpleSingleRegionControl.H"
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
    // Create the default single region mesh
    #include "createMesh.H"
    conductingRegionSolver regionSolver(runTime, mesh);
    solver& solver = regionSolver();

    // Create the outer PIMPLE loop and control structure
    pimpleSingleRegionControl pimple(solver.pimple);
    #include "createFields.H"

    // Set the initial time-step
    setDeltaT(runTime, solver);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "*** Preparing Elmer communications for sending" << nl << endl;
    Elmer<fvMesh> sending(mesh, //mesh
        1, // 1=send, -1=receive
        0, // 1=initialize, 0=w/o init
        1 // 1=multiregion/no O2E files, 0=exports O2E files
    );

    Info<< "*** Preparing Elmer communications for receiving" << nl << endl;
    Elmer<fvMesh> receiving(mesh, //mesh
        -1, // 1=send, -1=receive
        0, // 1=initialize, 0=w/o init
        1 // 1=multiregion/no O2E files, 0=exports O2E files
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
    //elmerClock = runTime.clockTimeIncrement();

    Info<< "\nStarting time loop\n" << endl;
    while (pimple.run(runTime) || lastTimeStep)
    {
        solver.preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solver);
        // Adjust time step so that last step is at end time.
        if (runTime.userTimeValue() + runTime.deltaTValue() > runTime.endTime().value())
        {
            const scalar lastDeltaT = runTime.endTime().value() - runTime.userTimeValue();
            runTime.setDeltaT(lastDeltaT);
            Info<< "Adjusting time step to match end time." << nl << endl;
            Info<< "deltaT = " << runTime.deltaTValue()  << nl << endl;
            lastTimeStep = true;
        }
        fieldPaths = getFieldPaths(mesh);
        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // PIMPLE corrector loop
        while (pimple.loop())
        {
            solver.moveMesh();
            solver.fvModels().correct();
            solver.prePredictor();
            solver.momentumPredictor();
            solver.thermophysicalPredictor();
            solver.pressureCorrector();
            solver.postCorrector();
        }

        solver.postSolve();
        //Update liquid-solid phase fraction
        /*if (solidificationEnabled)
        {
		    alpha1 = mesh().lookupObject<volScalarField>(solverSolidificationName);
        }*/
        //U = regionSolver.getVelocity();
        //p = regionSolver.getPressure();

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Check whether we need to update electromagnetic stuff with Elmer
        // Update electromagnetics with Elmer if magnetic field is
        // significantly disturbed and needs to be updated.
        // Otherwise, update electromagnetics locally by calculating electric potential.
        bool doElmer = regionSolver.updateMagneticField();

        // Calculate electric potential if current density will not be updated
        if (!doElmer)
        {
            regionSolver.setCorrectElectromagnetics();
        }
        regionSolver.solveElectromagnetics();
        const volVectorField& JxBRegion = regionSolver.getElectro().JxB;
        JxB = JxBRegion;
        const volScalarField& JJsigmaRegion = regionSolver.getElectro().JJsigma;
        JJsigma = JJsigmaRegion;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // Write last time step even if not write time
        if(runTime.writeTime() || lastTimeStep)
        {
            runTime.writeNow();
            // Cleanup
            forAll(fieldPaths, i)
            {
                if (!keepField[fieldPaths[i].first()])
                {
                    //Pout << "Deleting file " << fieldPaths[i].first() <<" was " <<
                    fileHandler().rm(fieldPaths[i].second());//<< endl;
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
