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

#include "fvCFD.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermophysicalTransportModel.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    //#include "createMesh.H"
    regionProperties rp(runTime);
    const wordList fluidNames
    (
        rp.found("fluid") ? rp["fluid"] : wordList(0)
    );
	
    // Populate fluid field pointer lists
    forAll(fluidNames, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidNames[i] << nl << endl;
        fvMesh fluidMesh
        (
            IOobject
            (
                fluidNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        );

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
        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            fluidMesh
        );
        
        T = (fluidMesh.C() & temperature_multiplier) +  temperature_addition;
        T.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
