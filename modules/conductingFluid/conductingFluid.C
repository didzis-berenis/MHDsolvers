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

#include "conductingFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(conductingFluid, 0);
    addToRunTimeSelectionTable(solver, conductingFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::conductingFluid::conductingFluid(fvMesh& mesh)
:
    isothermalFluid(mesh),

    electroBase(mesh),

    thermophysicalTransport
    (
        fluidThermoThermophysicalTransportModel::New
        (
            momentumTransport(),
            thermo
        )
    ),

    U_old_
    (
        IOobject
        (
            "U_old",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_
    ),

    U_old(U_old_)
{
    thermo.validate(type(), "h", "e");

    label PotERefCell = 0;
    scalar PotERefValue = 0.0;
    setRefCell
    ( 
        electro.PotE(),
        pimple.dict(),
        PotERefCell,
        PotERefValue
    );

    if (electro.isComplex())
    {
        setRefCell
        ( 
            electro.PotE(true),
            pimple.dict(),
            PotERefCell,
            PotERefValue
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::conductingFluid::~conductingFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::conductingFluid::prePredictor()
{
    isothermalFluid::prePredictor();

    if (pimple.predictTransport())
    {
        thermophysicalTransport->predict();
    }
}


void Foam::solvers::conductingFluid::postCorrector()
{
    isothermalFluid::postCorrector();

    if (pimple.correctTransport())
    {
        thermophysicalTransport->correct();
    }
}


void Foam::solvers::conductingFluid::solveElectromagnetics()
{
    if (electro.correctElectromagnetics())
    {
        // Update deltaU
        electroPtr->updateDeltaU(U_ - U_old_);
        //Correct current density
        electroPtr->correct();
        //Store old velocity for next update
        storeU();
    }
}


void Foam::solvers::conductingFluid::storeU()
{
    U_old_ = U_;
}

//non-const access for initialization purposes
Foam::volScalarField& Foam::solvers::conductingFluid::getTemperature()
{   
    return thermo_.T();
}

// ************************************************************************* //
