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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
	
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
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
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);
	/*
	forAll (T, cellI)
	{
		vector coords = mesh.C()[cellI];
		T[cellI] = (coords & temperature_multiplier).value() + temperature_addition.value();
	}
	*/
	T = (mesh.C() & temperature_multiplier) +  temperature_addition;
	T.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
