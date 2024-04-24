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
    checkRegionCellZones is based on OpenFOAM solvers and utilities

Description
    Checks if Pstream master thread has zero cells for any fluid region.

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include <set>
#include <fstream>
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    std::set<word> fluidNames = {"fluid", "conductingFluid", "incompressibleConductingFluid"};
    
    const fvMesh meshGlobal
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
    

    if (runTime.controlDict().found("regionSolvers"))
    {
        const dictionary& conductingRegionSolversDict =
            runTime.controlDict().subDict("regionSolvers");

        forAllConstIter(dictionary, conductingRegionSolversDict, iter)
        {
            const word regionName(iter().keyword());
            const word solverName(iter().stream());

            forAll(meshGlobal.cellZones(), zoneI)
            {
                const cellZone& cZone = meshGlobal.cellZones()[zoneI];

                bool hasZeroCells = cZone.name() == regionName
                && fluidNames.find(solverName) != fluidNames.end()
                && cZone.size() == 0;

                if (hasZeroCells)
                {
                    std::ofstream outFile("HAS_ZERO_CELLS", std::ofstream::out);
                    outFile.close();
                }

                bool masterHasZeroCells = Pstream::master()
                && hasZeroCells;

                if (masterHasZeroCells)
                {
                    std::ofstream outFile("MASTER_HAS_ZERO_CELLS", std::ofstream::out);
                    outFile.close();
                }
            }

        }
    }
    else
    {
        FatalIOErrorInFunction(runTime.controlDict())
            << "regionSolvers list missing from "
            << runTime.controlDict().name()
            << exit(FatalIOError);
    }

    return 0;
}


// ************************************************************************* //
