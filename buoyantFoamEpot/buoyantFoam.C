/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Original buoyantFoam solver is part of OpenFOAM.

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
    buoyantFoamMHD is a modified buoyantFoam solver, where the modification 
    is based on EOF-Library solver mhdVxBPimpleFoam.

Description
    Solver for steady or transient buoyant, turbulent flow of compressible
    fluids for electromagnetically forced and heated flows, with optional 
    mesh motion and mesh topology changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermophysicalTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "hydrostaticInitialisation.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "Elmer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	const Foam::vector ivec = Foam::vector(1,0,0);
	const Foam::vector jvec = Foam::vector(0,1,0);
	const Foam::vector kvec = Foam::vector(0,0,1);
    
	const pointField& cellCentres = mesh.cellCentres();
	
	double xDim = std::abs(max(cellCentres & ivec ) - min(cellCentres & ivec ));
	double yDim = std::abs(max(cellCentres & jvec ) - min(cellCentres & jvec ));
	double zDim = std::abs(max(cellCentres & kvec ) - min(cellCentres & kvec ));
	
	double Lchar = 0.5*std::min(xDim,std::min(yDim,zDim));
	
	double Rem0 = 4*3.14159*(std::pow(10,-7))*sigma.value()*Lchar;
    Info<< "Rem0 = " << Rem0 << endl;

    double OFClock = 0;
    double elmerClock = runTime.clockTimeIncrement();

    Info<< "\nStarting time loop\n" << endl;

    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    sending.sendVector(U);

    // Receive fields from Elmer
    Elmer<fvMesh> receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(Jre);
    receiving.recvVector(Jim);
    receiving.recvVector(Bre);
    receiving.recvVector(Bim);
    
	//Lorentz force term initialization
	JxB =  0.5*((Jre ^ Bre) + (Jim ^ Bim) );
	JJsigma =  0.5*((Jre . Jre) + (Jim . Jim) )/sigma;
	
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

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        fvModels.preUpdateMesh();

        // Store momentum to set rhoUf for introduced faces.
        autoPtr<volVectorField> rhoU;
        if (rhoUf.valid())
        {
            rhoU = new volVectorField("rhoU", rho*U);
        }

        // Update the mesh for topology change, mesh to mesh mapping
        mesh.update();

        runTime++;

        Info<< "SimulationTime = " << runTime.userTimeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (!pimple.flow())
            {
                if (pimple.models())
                {
                    fvModels.correct();
                }

                if (pimple.thermophysics())
                {
                    #include "EEqn.H"
                }
            }
            else
            {
                if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
                {
                    // Move the mesh
                    mesh.move();

                    if (mesh.changing())
                    {
                        gh = (g & mesh.C()) - ghRef;
                        ghf = (g & mesh.Cf()) - ghRef;

                        MRF.update();

                        if (correctPhi)
                        {
                            #include "correctPhi.H"
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

                if (pimple.firstPimpleIter() && !pimple.simpleRho())
                {
                    #include "rhoEqn.H"
                }

                if (pimple.models())
                {
                    fvModels.correct();
                }

                #include "UEqn.H"

                if (pimple.thermophysics())
                {
                    #include "EEqn.H"
                }

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                    thermophysicalTransport->correct();
                }
            }
        }

        rho = thermo.rho();
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Check whether we need to update electromagnetic stuff with Elmer

        dimensionedScalar smallU
        (
            "smallU",
            dimensionSet(0, 1, -1, 0, 0, 0 ,0),
            1e-6
        );

        bool doElmer = false;
        
        scalar maxRemDiff_local = Rem0*max(mag(U_old-U)).value();        
        
        scalar maxRelDiff_local = (max(mag(U_old-U)/(average(mag(U))+smallU))).value();
        
        if((maxRelDiff_local>maxRelDiff || maxRelDiff<SMALL) && maxRelDiff+SMALL<=1.0) {
            doElmer = true;
        }
        else if(maxRemDiff_local>maxRemDiff && maxRelDiff-SMALL<=1.0) {
            doElmer = true;
        }

        // Calculate electric potential if current density will not be updated
        if (!doElmer)
        {
            volVectorField JUBre = Jre;
            {
                #include "PotEreEqn.H"
            }
            volVectorField JUBim = Jim;
            {
                #include "PotEimEqn.H"
            }
            JxB =  0.5*(((Jre+JUBre) ^ Bre) + ((Jim+JUBim) ^ Bim) );
	        JJsigma =  0.5*(((Jre+JUBre) . (Jre+JUBre)) + ((Jim+JUBim) . (Jim+JUBim)) )/sigma;
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        runTime.write();
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
            receiving.recvVector(Jre);
            receiving.recvVector(Jim);
            receiving.recvVector(Bre);
            receiving.recvVector(Bim);
            JxB =  0.5*((Jre ^ Bre) + (Jim ^ Bim) );
	        JJsigma =  0.5*((Jre . Jre) + (Jim . Jim) )/sigma;
			
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
    receiving.recvVector(Jre);
    receiving.recvVector(Jim);
    receiving.recvVector(Bre);
    receiving.recvVector(Bim);
	
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
