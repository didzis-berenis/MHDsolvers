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

    Info<< "\nStarting time loop\n" << endl;

    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    sending.sendVector(U);

    // Receive fields from Elmer
    Elmer<fvMesh> receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(JxB);
    receiving.recvScalar(JJsigma);
/*
	//alternatively receive current and magnetic field separately
    receiving.recvVector(J);
    receiving.recvVector(B);
    JxB = J ^ B;
    JJsigma = (J & J)/sigma;
*/

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

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

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

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
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
			receiving.recvScalar(JJsigma);
			/*
			//alternatively receive current and magnetic field separately
    		receiving.recvVector(J);
    		receiving.recvVector(B);
    		JxB = J ^ B;
    		JJsigma = (J & J)/sigma;
			*/
        }
    }

    Info << "Final iteration: "
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

    //Final iter for Elmer
    U_old = U;
    // Send fields to Elmer
    sending.sendStatus(0);
    sending.sendVector(U);
    // Receive fields form Elmer
    receiving.sendStatus(0);
   	receiving.recvVector(JxB);
    receiving.recvScalar(JJsigma);
	/*
	//alternatively receive current and magnetic field separately
    receiving.recvVector(J);
    receiving.recvVector(B);
    JxB = J ^ B;
    JJsigma = (J & J)/sigma;
	*/

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
