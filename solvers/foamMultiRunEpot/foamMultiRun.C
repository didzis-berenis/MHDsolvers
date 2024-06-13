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
    modification was made to update electrical currents in OpenFOAM, while the 
    change in magnetic Reynolds number doesn't exceed the provided value. This 
	modification was based on the epotFoam solver, which can be found
	in https://doi.org/10.13140/RG.2.2.12839.55201 (Chapter 4).

Description
    Solver for steady or transient electromagnetically forced fluid flow and 
    solid heat conduction, with conjugate heat transfer between regions, 
    buoyancy effects, turbulence and radiation modelling.

    Compile option ELMER_TIME == HARMONIC_TIME builds foamMultiRunEpot
    solver, which assumes coupling with harmonic (time-averaged) ElmerFEM solver.

    Compile option ELMER_TIME == TRANSIENT_TIME builds foamMultiRunEpotTransient
    solver, which assumes coupling with transient ElmerFEM solver.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "conductingRegionSolvers.H"
#include "pimpleMultiRegionControl.H"
#include "setDeltaT.H"
#include "findRefCell.H"

using namespace Foam;
#include "Elmer.H"
#include <fstream>
#include "globalRegionMapper.H"
#include "fieldMapper.H"
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

    // Create the region meshes and solvers
    conductingRegionSolvers solvers(runTime);
    #include "createMeshes.H"
    // Create the outer PIMPLE loop and control structure
    pimpleMultiRegionControl pimple(runTime, solvers);
    #include "createFields.H"

    // Set the initial time-step
    setDeltaT(runTime, solvers);
    #if (ELMER_TIME == HARMONIC_TIME)
        #include "createHarmonicEpotControls.H"
    #elif (ELMER_TIME == TRANSIENT_TIME)
        #include "createTransientEpotControls.H"
    #endif
    #include "createElmerComms.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    double OFClock = 0;
    double elmerClock = runTime.clockTimeIncrement();

    Info<< nl << "Starting time loop\n" << endl;

    bool initialize_elmer = true;
    int elmer_status = 1; // 1=ok, 0=lastIter, -1=error
    #if (ELMER_TIME == HARMONIC_TIME)
        #include "setHarmonicElmerComms.H"
    #elif (ELMER_TIME == TRANSIENT_TIME)
        #include "setTransientElmerComms.H"
    #endif
    initialize_elmer = false;
	
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

    // Write initial values
    #include "writeIntegrals.H"

    while (pimple.run(runTime))
    {
        forAll(solvers, i)
        {
            solvers[i].preSolve();
        }

        solvers.setGlobalPrefix();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(runTime, solvers);
        // Update paths for cleanup
        forAll(regionNames, i)
        {
            regionPaths[i] = getRegionPath(solvers.mesh(regionNames[i]));
        }
        fieldPaths = getFieldPaths(meshGlobal);
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

        solvers.setGlobalPrefix();
        //Update liquid-solid phase fraction
        forAll(regionNames, i)
        {
            if (solidificationEnabled[i])
            {
                alpha1Region[i] = solvers.mesh(regionNames[i]).lookupObject<volScalarField>(solverSolidificationName);
                scalarFieldToGlobal(alpha1Global,alpha1Region[i],regionNames[i]);
            }
        }
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
        #include "writeGlobal.H"
        #include "writeIntegrals.H"
        OFClock = runTime.clockTimeIncrement();

        Info<< "ExecutionTime : " << "Hydrodynamics step = " << OFClock << " s"
            << " ; Electrodynamics step = " << elmerClock << " s"
            << " ; ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
			
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Update electromagnetic stuff with Elmer

        if(doElmer && runTime.run()) {
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
    elmer_status = 0; // 1=ok, 0=lastIter, -1=error
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
