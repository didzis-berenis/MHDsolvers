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
//#include "fluidThermo.H"
//#include "compressibleMomentumTransportModels.H"
//#include "fluidThermophysicalTransportModel.H"

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
        //const word& solverName = solverNames[i].second();
        if (solvers.isFluid(regionName))
        {
            Info<< "*** Reading fluid properties for region "
                << regionName << nl << endl;

            fvMesh& fluidMesh = solvers.mesh(regionName);

            IOdictionary physicalProperties
                (
                    IOobject
                    (
                        "physicalProperties",
                        runTime.constant(),
                        fluidMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

            dimensionedVector temperature_multiplier
            (
                "temperature_multiplier",
                dimensionSet(0, -1, 0, 1, 0, 0, 0),
                physicalProperties
            );
            dimensionedScalar temperature_addition
            (
                "temperature_addition",
                dimensionSet(0, 0, 0, 1, 0, 0, 0),
                physicalProperties
            );
            
            Info<< "Reading field T\n" << endl;
            volScalarField& T = solvers.getTemperature(regionName);
            
            T = (fluidMesh.C() & temperature_multiplier) +  temperature_addition;
            T.write();

        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
