/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    chtMultiRegionFoamEpotTransient is a modified chtMultiRegionFoam solver, 
    where the modification is based on EOF-Library solver mhdVxBPimpleFoam. 
    Additional modification was made to update electrical currents in OpenFOAM, 
    while the change in magnetic Reynolds number doesn't exceed the provided 
    value. This modification was based on the epotFoam solver, which can be 
    found in https://doi.org/10.13140/RG.2.2.12839.55201 (Chapter 4).

Description
    Solver for steady or transient electromagnetically forced fluid flow and 
    solid heat conduction, with conjugate heat transfer between regions, 
    buoyancy effects, turbulence, reactions and radiation modelling. 
    buoyantFoamEpot assumes coupling with transient ElmerFEM solver.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "pressureReference.H"
#include "hydrostaticInitialisation.H"
#include "Elmer.H"
#include "globalRegionInterpolator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createFluidPressureControls.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"
    #include "createElmerComms.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	int writeMultiplier( readScalar(runTime.controlDict().lookup("writeMultiplier")) );
	long long int writeCounter = 0;
	
	word writeControlDict(runTime.controlDict().lookup("writeControl"));
	const bool adjustableRunTime = (writeControlDict=="adjustableRunTime");

    double OFClock = 0;
    double elmerClock = runTime.clockTimeIncrement();

    Info<< "\nStarting time loop\n" << endl;

    bool initialize_elmer = true;
    #include "setElmerComms.H"
    initialize_elmer = false;
	
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

    while (pimples.run(runTime))
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Optional number of energy correctors
        const int nEcorr = pimples.dict().lookupOrDefault<int>
        (
            "nEcorrectors",
            1
        );

        // --- PIMPLE loop
        while (pimples.loop())
        {
            List<tmp<fvVectorMatrix>> UEqns(fluidRegions.size());

            for(int Ecorr=0; Ecorr<nEcorr; Ecorr++)
            {
                forAll(solidRegions, i)
                {
                    Info<< "\nSolving for solid region "
                        << solidRegions[i].name() << endl;
                    #include "setRegionSolidFields.H"
                    #include "solveSolid.H"
                }
                forAll(fluidRegions, i)
                {
                    Info<< "\nSolving for fluid region "
                        << fluidRegions[i].name() << endl;
                    #include "setRegionFluidFields.H"
                    #include "solveFluid.H"
                }
            }
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Check whether we need to update electromagnetic stuff with Elmer

        bool doElmer = false;
        
        if (adjustableRunTime)
		{
			if (runTime.writeTime()) doElmer = true;
		}
		else
		{
			writeCounter++;
			if ( (writeCounter % writeMultiplier) == 0 && runTime.run()) doElmer = true;
		}

        // Calculate electric potential if current density will not be updated
        if (!doElmer)
        {
            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];
                volVectorField& U = UFluid[i];
                volVectorField& U_old = U_oldFluid[i];
                volVectorField& J = JFluid[i];
                volVectorField& B = BFluid[i];
                dimensionedScalar& sigma = sigmaFluid[i];
                volScalarField& PotE = PotEFluid[i];

                volVectorField JUB = J;
                {
                    #include "PotEEqn.H"
                }
                JxBFluid[i] =  ((J+JUB) ^ B );
                JJsigmaFluid[i] =  ((J+JUB) & (J+JUB))/sigma;
            }
        }

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
            #include "setElmerComms.H"
			
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
    #include "setElmerComms.H"
	
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
