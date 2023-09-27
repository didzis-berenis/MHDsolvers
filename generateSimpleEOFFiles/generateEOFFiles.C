/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Original simpleFoam solver is part of OpenFOAM.

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
    generatePimpleEOFFiles is based on code in EOF-Library solver mhdVxBPimpleFoam.

Description
    Case is expected to be initialized for SIMPLE solver. 
	U field is sent to Elmer once to generate O2E files.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "Elmer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
	
	Info<< "Reading field U\n" << endl;
	volVectorField U
	(
		IOobject
		(
			"U",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);

    Info<< "\nSend fields to Elmer once\n" << endl;

    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(0); // 1=ok, 0=lastIter, -1=error
    sending.sendVector(U);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
