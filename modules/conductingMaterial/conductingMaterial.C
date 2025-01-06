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

\*---------------------------------------------------------------------------*/

#include "conductingMaterial.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(conductingMaterial, 0);
    addToRunTimeSelectionTable(solver, conductingMaterial, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::conductingMaterial::conductingMaterial
(
    fvMesh& mesh
)
:
    solver(mesh),
    electroBase(mesh)
{

/*label PotERefCell = 0;
scalar PotERefValue = 0.0;
setRefCell
( 
    electro.PotE(),
    pimple.dict(),
    PotERefCell,
    PotERefValue
);
mesh.schemes().setFluxRequired(electro.PotE().name());

if (electro.isComplex())
{
    setRefCell
    ( 
        electro.PotE(true),
        pimple.dict(),
        PotERefCell,
        PotERefValue
    );
    mesh.schemes().setFluxRequired(electro.PotE(true).name());
}*/

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::conductingMaterial::~conductingMaterial()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalar Foam::solvers::conductingMaterial::maxDeltaT() const
{
    scalar deltaT = fvModels().maxDeltaT();
    return deltaT;
}


void Foam::solvers::conductingMaterial::preSolve()
{}


void Foam::solvers::conductingMaterial::moveMesh()
{}


void Foam::solvers::conductingMaterial::prePredictor()
{}


void Foam::solvers::conductingMaterial::momentumPredictor()
{}


void Foam::solvers::conductingMaterial::thermophysicalPredictor()
{}


void Foam::solvers::conductingMaterial::pressureCorrector()
{}


void Foam::solvers::conductingMaterial::postCorrector()
{}


void Foam::solvers::conductingMaterial::postSolve()
{}


/*void Foam::solvers::conductingMaterial::postCorrector()
{
    if (pimple.correctTransport())
    {
        thermophysicalTransport->correct();
    }
    if (electro.correctElectromagnetics())
    {
        //Correct current density
        electro_.correct();
    }
}*/
//non-const access for initialization purposes
/*Foam::volScalarField& Foam::solvers::conductingMaterial::getTemperature()
{   
    return thermo_.T();
}*/

void Foam::solvers::conductingMaterial::solveElectromagnetics()
{
    electro_.solve();
}

// ************************************************************************* //
