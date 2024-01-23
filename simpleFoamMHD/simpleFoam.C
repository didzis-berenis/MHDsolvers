/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Original simpleFoam solver is part of OpenFOAM.

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
    simpleFoamMHD is a modified pimpleFoam solver, where the modification 
    is based on EOF-Library solver mhdVxBPimpleFoam.

Description
    Steady-state solver for electromagnetically forced incompressible, 
    turbulent flow, using the SIMPLE algorithm.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModels.H"
#include "simpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "Elmer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

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
    receiving.recvVector(JxB);
	
    // Create file for logging simulation times whenever Elmer is called
    string elmerTimesFileName = "elmerTimes.log";
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
    
    while (simple.loop(runTime))
    {
        Info<< "SimulationTime = " << runTime.userTimeName() << nl << endl;

        fvModels.correct();

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        viscosity->correct();
        turbulence->correct();

        runTime.write();
        OFClock = runTime.clockTimeIncrement();

        Info<< "ExecutionTime : " << "Hydrodynamics step = " << OFClock << " s"
            << " ; Electrodynamics step = " << elmerClock << " s"
            << " ; ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
			
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        dimensionedScalar smallU
        (
            "smallU",
            dimensionSet(0, 1, -1, 0, 0, 0 ,0),
            1e-6
        );

        // Check whether we need to update electromagnetic stuff with Elmer
        scalar maxRelDiff_local = (max(mag(U_old-U)/(average(mag(U))+smallU))).value();

        bool doElmer = false;
        if(maxRelDiff_local>maxRelDiff && (maxRelDiff<SMALL || maxRelDiff+SMALL<=1.0)) {
            doElmer = true;
        }

        if(doElmer && runTime.run()) {
            U_old = U;

            // Send fields to Elmer
            sending.sendStatus(1);
            sending.sendVector(U);

            // Receive fields form Elmer
            receiving.sendStatus(1);
    	    receiving.recvVector(JxB);
			
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
   	receiving.recvVector(JxB);
	
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
