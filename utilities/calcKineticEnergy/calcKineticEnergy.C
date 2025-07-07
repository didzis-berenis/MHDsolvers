/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    calcKineticEnergy is based on OpenFOAM solvers and utilities

Description
    Calculates kinetic energy for a decomposed case

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "processorRunTimes.H"
#include "domainDecomposition.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "reconstructLagrangian.H"

using namespace Foam;
#include <fstream>
#include "fvFieldReconstructor.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, true);
    argList::noParallel();
    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
    #include "setRootCase.H"

    // Set time from database
    Info<< "Create time\n" << endl;
    processorRunTimes runTimes(Foam::Time::controlDictName, args);

    // Allow override of time
    const instantList times = runTimes.selectProc(args);

    const Time& runTime = runTimes.procTimes()[0];

    #include "setRegionNames.H"

    // Determine the processor count
    const label nProcs = fileHandler().nProcs
    (
        args.path(),
        regionNames[0] == polyMesh::defaultRegion
      ? word::null
      : regionNames[0]
    );

    if (!nProcs)
    {
        FatalErrorInFunction
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Warn fileHandler of number of processors
    const_cast<fileOperation&>(fileHandler()).setNProcs(nProcs);

    if (times.empty())
    {
        WarningInFunction << "No times selected" << endl;
        exit(1);
    }

    // Reconstruct all regions
    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];//singleRegion = "region0"

        // Create meshes
        Info<< "\n\nReconstructing mesh " << regionName << nl << endl;
        domainDecomposition meshes(runTimes, regionName);
        meshes.readReconstruct(true);
        
        fvFieldReconstructor fvReconstructor
        (
            meshes.completeMesh(),
            meshes.procMeshes(),
            meshes.procFaceAddressing(),
            meshes.procCellAddressing(),
            meshes.procFaceAddressingBf()
        );

        Info<< "Starting time loop " << nl << endl;
        // Loop over all times
        forAll(times, timei)
        {
            // Set the time
            runTimes.setTime(times[timei], timei);

            Info<< "Time = " << runTimes.completeTime().userTimeName()
                << nl << endl;

            /*const word fieldName = "U";
            const IOobject object
            (
                fieldName,
                runTimes.completeTime().userTimeName(),
                meshes.procMeshes()[0]
            );*/
            const volVectorField U(fvReconstructor.reconstructFvVolumeField<vector>(
                IOobject
                (
                    word("U"),
                    runTimes.completeTime().userTimeName(),
                    meshes.procMeshes()[0]
                )
            )());
            #include "writeIntegrals.H"
        }
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
