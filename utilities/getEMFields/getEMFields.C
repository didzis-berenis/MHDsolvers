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

using namespace Foam;
#include "Elmer.H"
#include <fstream>
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

    // Initialize fields
    //volVectorField JxB = regionSolver.getElectro().JxB;
    //volScalarField JJsigma = regionSolver.getElectro().JJsigma;
    volVectorField Jre = regionSolver.getElectro().J();
    volVectorField Jim = regionSolver.getElectro().J(true);
    volVectorField Bre = regionSolver.getElectro().B();
    volVectorField Bim = regionSolver.getElectro().B(true);
    const volVectorField& U = regionSolver.getU();
    const volScalarField& sigma = regionSolver.getElectro().sigma();
    // Cannot use const reference for Elmer functions,
    // so preparing designated field for sending velocity to Elmer.
    volVectorField Usent = U;
    volScalarField sigmaSent = sigma;

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

    Info<< nl << "Initializing electromagnetic solver\n" << endl;

    int elmer_status = 1; // 1=ok, 0=lastIter, -1=error
    bool initialize_elmer = true;
    int nElmerCorrectors = (runTime.controlDict().lookupOrDefault("nElmerCorrectors",1));
    Info<< "Initializing electromagnetic solver" << nl << endl;
    //Final iter for Elmer
    elmer_status = 0; // 1=ok, 0=lastIter, -1=error

    // Save old U field
    regionSolver.storeU();

    // Save old phase-fraction field
    regionSolver.storeAlpha1();

    for (int iter = 0; iter < nElmerCorrectors; iter++)
    {
        volVectorField pJre = Jre;
        volVectorField pJim = Jim;
        volVectorField pBre = Bre;
        volVectorField pBim = Bim;
        Info << "\nElmer iteration: " << iter << endl;
        // Send fields to Elmer
        Info << "Sending fields to Elmer" << endl;
        if (initialize_elmer) sending.initialize();
        sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
        //cannot use const reference, so assign to new field 
        Usent = U;
        sigmaSent = sigma;
        sending.sendVector(Usent);
        //if (regionSolver.isIncompressibleConductingVoF())
        {
            sending.sendScalar(sigmaSent);
        }
        // Receive fields from Elmer
        Info << "Receiving fields from Elmer" << endl;
        if (initialize_elmer) receiving.initialize();
        receiving.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
        receiving.recvVector(Jre);
        if (regionSolver.isElectroHarmonic())
        {
            receiving.recvVector(Jim);
        }
        receiving.recvVector(Bre);
        if (regionSolver.isElectroHarmonic())
        {
            receiving.recvVector(Bim);
        }
        scalar eJre = gMax(mag(pJre - Jre)())/(gMax(mag(pJre)())+SMALL);
        scalar eJim = gMax(mag(pJim - Jim)())/(gMax(mag(pJim)())+SMALL);
        scalar eBre = gMax(mag(pBre - Bre)())/(gMax(mag(pBre)())+SMALL);
        scalar eBim = gMax(mag(pBim - Bim)())/(gMax(mag(pBim)())+SMALL);
        scalar maxError = max(max(eJre,eJim),max(eBre,eBim));
        bool quitLoop = maxError < ROOTSMALL;
        // Initialization is needed only on the first update.
        if (initialize_elmer)
            initialize_elmer = false;
        //Quit the loop if bad status
        quitLoop |= elmer_status != 1;
        if (quitLoop) break;
    }
    //update solver fields
    regionSolver.setJ(Jre);
    regionSolver.setB(Bre);
    if (regionSolver.isElectroHarmonic())
    {
        regionSolver.setJ(Jim,true);
        regionSolver.setB(Bim,true);
    }

    runTime.writeNow();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
