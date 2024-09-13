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
    initializeTemperature is based on OpenFOAM solvers and utilities

Description
    Case is expected to be initialized with temperature distribution

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "conductingRegionSolvers.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // Create the region meshes and solvers
    conductingRegionSolvers solvers(runTime);
    List<Pair<word>> solverNames = solvers.getNames();

    forAll(solverNames, i)
    {
        const word& regionName = solverNames[i].first();
        solvers.calcTemperatureGradient(regionName);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
