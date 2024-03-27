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
    foamRun

Description
    Loads and executes an OpenFOAM solver module either specified by the
    optional \c solver entry in the \c controlDict or as a command-line
    argument.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

Usage
    \b foamRun [OPTION]

      - \par -solver <name>
        Solver name

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

    Example usage:
      - To run a \c rhoPimpleFoam case by specifying the solver on the
        command line:
        \verbatim
            foamRun -solver fluid
        \endverbatim

      - To update and run a \c rhoPimpleFoam case add the following entries to
        the controlDict:
        \verbatim
            application     foamRun;

            solver          fluid;
        \endverbatim
        then execute \c foamRun

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "conductingRegionSolver.H"
#include "pimpleSingleRegionControl.H"
#include "setDeltaT.H"
#include "findRefCell.H"

using namespace Foam;
#include "Elmer.H"
#include <fstream>
#define TRANSIENT_TIME  2
#define HARMONIC_TIME   3
#if (ELMER_TIME == HARMONIC_TIME)
#warning "Compiling for coupling with HARMONIC Elmer simulation!"
#elif (ELMER_TIME == TRANSIENT_TIME)
#warning "Compiling for coupling with TRANSIENT Elmer simulation!"
#else
#error "Please define appropriate functions for your Elmer simulation!"
#endif

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

    double OFClock = 0;
    double elmerClock = runTime.clockTimeIncrement();

    Info<< "\nStarting time loop\n" << endl;

    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh, //mesh
        1, // 1=send, -1=receive
        1, // 1=initialize, 0=w/o init
        1 // 1=multiregion/no O2E files, 0=exports O2E files
    );
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    sending.sendVector(U);

    // Receive fields from Elmer
    Elmer<fvMesh> receiving(mesh, //mesh
        -1, // 1=send, -1=receive
        1, // 1=initialize, 0=w/o init
        1 // 1=multiregion/no O2E files, 0=exports O2E files
    );
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    #if (ELMER_TIME == HARMONIC_TIME)
        #include "setHarmonicElmerComms.H"
    #elif (ELMER_TIME == TRANSIENT_TIME)
        #include "setTransientElmerComms.H"
    #endif
	
    // Create file for logging simulation times whenever Elmer is called
    string elmerTimesFileName = "postProcessing/elmerTimes.log";
	// Log the current simulation time
    if (Pstream::master())
    {
		std::ofstream elmerTimes(elmerTimesFileName, std::ios::app);
		if (elmerTimes.is_open())
		{
			elmerTimes << runTime.timeName() << std::endl;
			elmerTimes.close();
		}
		else FatalErrorInFunction << "ERROR: Couldn't open " << elmerTimesFileName << " for writing!\n" << abort(FatalError);
	}

    elmerClock = runTime.clockTimeIncrement();

    while (pimple.run(runTime))
    {
        solver.preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solver);

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
        #include "getSolverFields.H"
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Check whether we need to update electromagnetic stuff with Elmer
        bool doElmer = false;

        #if (ELMER_TIME == HARMONIC_TIME)
            #include "setHarmonicPotential.H"
        #elif (ELMER_TIME == TRANSIENT_TIME)
            #include "setTransientPotential.H"
        #endif

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        runTime.write();
        #include "writeIntegrals.H"
        OFClock = runTime.clockTimeIncrement();

        Info<< "ExecutionTime : " << "Hydrodynamics step = " << OFClock << " s"
            << " ; Electrodynamics step = " << elmerClock << " s"
            << " ; ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
			
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update electromagnetic stuff with Elmer

        if(doElmer && runTime.run()) {
            U_old = U;

            // Send fields to Elmer
            sending.sendStatus(1);
            sending.sendVector(U);

            // Receive fields form Elmer
            receiving.sendStatus(1);
            #if (ELMER_TIME == HARMONIC_TIME)
                #include "setHarmonicElmerComms.H"
            #elif (ELMER_TIME == TRANSIENT_TIME)
                #include "setTransientElmerComms.H"
            #endif
			
			// Log the current simulation time
			if (Pstream::master())
			{
				std::ofstream elmerTimes(elmerTimesFileName, std::ios::app);
				if (elmerTimes.is_open())
				{
					elmerTimes << runTime.timeName() << std::endl;
					elmerTimes.close();
				}
				else FatalErrorInFunction << "ERROR: Couldn't open " << elmerTimesFileName << " for writing!\n" << abort(FatalError);
			}
        }
        elmerClock = runTime.clockTimeIncrement();
    }

    Info << "Final iteration: "
        << " ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

    //Final iter for Elmer
    U_old = U;
    // Send fields to Elmer
    sending.sendStatus(0);
    sending.sendVector(U);
    // Receive fields form Elmer
    receiving.sendStatus(0);
    #if (ELMER_TIME == HARMONIC_TIME)
        #include "setHarmonicElmerComms.H"
    #elif (ELMER_TIME == TRANSIENT_TIME)
        #include "setTransientElmerComms.H"
    #endif
	
	// Log the current simulation time
    if (Pstream::master())
    {
		std::ofstream elmerTimes(elmerTimesFileName, std::ios::app);
		if (elmerTimes.is_open())
		{
			elmerTimes << runTime.timeName() << std::endl;
			elmerTimes.close();
		}
		else FatalErrorInFunction << "ERROR: Couldn't open " << elmerTimesFileName << " for writing!\n" << abort(FatalError);
	}

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
