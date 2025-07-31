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
#include "fvcCurl.H"

using namespace Foam;
#include "Elmer.H"
#include <fstream>
#include "fieldMapper.H"
const scalar ALMOST_ONE = 1.0 - 1e-6;
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

    int elmer_status = 0; // 1=ok, 0=lastIter, -1=error
    bool initialize_elmer = true;
    Info<< "Initializing electromagnetic solver" << nl << endl;
    #include "runElmerUpdate.H"
    scalarField volume = mesh.V();
    scalar totalVolume = gSum(volume);
    volVectorField JxB = regionSolver.getElectro().JxB;
    volScalarField JxB_X((JxB & Foam::vector(1,0,0)));
    volScalarField JxB_Y((JxB & Foam::vector(0,1,0)));
    volScalarField JxB_Z((JxB & Foam::vector(0,0,1)));
    scalar JxB_X_Integral = gSum(JxB_X*volume)/totalVolume;
    scalar JxB_Y_Integral = gSum(JxB_Y*volume)/totalVolume;
    scalar JxB_Z_Integral = gSum(JxB_Z*volume)/totalVolume;
    scalar JxB_X_IntegralAbs = gSum(mag(JxB_X)*volume)/totalVolume;
    scalar JxB_Y_IntegralAbs = gSum(mag(JxB_Y)*volume)/totalVolume;
    scalar JxB_Z_IntegralAbs = gSum(mag(JxB_Z)*volume)/totalVolume;

    Info << "JxB integral = ( " << JxB_X_Integral << "\t" << JxB_Y_Integral << "\t" << JxB_Z_Integral << " )" << endl;
    Info << "JxB integral mag = ( " << JxB_X_IntegralAbs << "\t" << JxB_Y_IntegralAbs << "\t" << JxB_Z_IntegralAbs << " )" << endl;

    volVectorField curlU(fvc::curl(JxB));
    volScalarField curlUx((curlU & Foam::vector(1,0,0)));
    volScalarField curlUy((curlU & Foam::vector(0,1,0)));
    volScalarField curlUz((curlU & Foam::vector(0,0,1)));
    scalar curlUxIntegral = gSum(curlUx*volume)/totalVolume;
    scalar curlUyIntegral = gSum(curlUy*volume)/totalVolume;
    scalar curlUzIntegral = gSum(curlUz*volume)/totalVolume;
    scalar curlUxIntegralAbs = gSum(mag(curlUx)*volume)/totalVolume;
    scalar curlUyIntegralAbs = gSum(mag(curlUy)*volume)/totalVolume;
    scalar curlUzIntegralAbs = gSum(mag(curlUz)*volume)/totalVolume;

    Info << "curl(JxB) integral = ( " << curlUxIntegral << "\t" << curlUyIntegral << "\t" << curlUzIntegral << " )" << endl;
    Info << "curl(JxB) integral mag = ( " << curlUxIntegralAbs << "\t" << curlUyIntegralAbs << "\t" << curlUzIntegralAbs << " )" << endl;

    runTime.writeNow();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
