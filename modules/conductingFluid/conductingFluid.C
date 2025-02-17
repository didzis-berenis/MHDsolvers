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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::conductingFluid::readControls()
{
    isothermalFluid::readControls();

    maxDi =
        runTime.controlDict().lookupOrDefault<scalar>("maxDi", 1.0);
}


void Foam::solvers::conductingFluid::correctDiNum()
{
    const volScalarField kappa = thermo.kappa();
    const scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()/rho.primitiveField()
    );
    // Alpha/rho from material properties
    const scalarField sumAlpha
    (
        fvc::surfaceSum
        (
            mesh.magSf()
        *fvc::interpolate(kappa)
        *mesh.surfaceInterpolation::deltaCoeffs()
        )()()/(thermo.rho()()*thermo.Cp()())
    );
    // Alpha/rho from turbulence
    const scalarField sumAlphat
    (
        alphat.primitiveField()/rho.primitiveField()
    );
    // Diffusion number including turbulence
    const scalarField DiNumvf
    (
        (sumAlpha+sumAlphat)*runTime.deltaT().value()/mesh.V()
    );

    const scalar meanDiNum = gAverage(DiNumvf);
    const scalar maxDiNum = gMax(DiNumvf);

    Info<< "Diffusion Number mean: " << meanDiNum
        << " max: " << maxDiNum << endl;

    DiNum = maxDiNum;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::conductingFluid::conductingFluid(fvMesh& mesh)
:
    isothermalFluid(mesh),

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

    U_old(U_old_),

    electroBase(mesh),

    alphat_
    (
            IOobject
            (
            "alphat",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimVelocity*dimDensity*dimLength,0)
    ),

    alphat
    (
        mesh.objectRegistry::foundObject<volScalarField>("alphat") ?
        mesh.objectRegistry::lookupObjectRef<volScalarField>("alphat") :
        alphat_
    )
{
    thermo.validate(type(), "h", "e");
    if (transient())
    {
        correctDiNum();
    }

    label PotERefCell = 0;
    scalar PotERefValue = 0.0;
    setRefCell
    ( 
        electro.PotE(),
        pimple.dict(),
        PotERefCell,
        PotERefValue
    );

    if (electro_.isComplex())
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

Foam::scalar Foam::solvers::conductingFluid::maxDeltaT() const
{
    // Combine hydrodynamic and thermophysical conditions for maxDeltaT
    scalar deltaT = isothermalFluid::maxDeltaT();

    if (DiNum > small)
    {
        deltaT = min(deltaT, maxDi/DiNum*runTime.deltaTValue());
    }

    return deltaT;
}


void Foam::solvers::conductingFluid::preSolve()
{
    // Combine hydrodynamic and thermophysical preSolve actions
    isothermalFluid::preSolve();

    // Read the controls
    readControls();

    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();

    if (transient())
    {
        correctDiNum();
    }
}


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
    if (electro.correctElectromagnetics())
    {
        //Correct current density
        electro_.correct();
    }
}


void Foam::solvers::conductingFluid::solveElectromagnetics()
{
    //Solve potential equation
    electro_.solve();
}


void Foam::solvers::conductingFluid::storeU()
{
    U_old_ = U_;
    //Store U_old_ boundary field values if needed
    /*volVectorField::Boundary& U_oldBf = U_old_.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U_.boundaryField();
    forAll(U_oldBf, patchi)
    {
        fvPatchVectorField& pU_old = U_oldBf[patchi];
        pU_old = UBf[patchi];
    }*/
}

//non-const access for initialization purposes
Foam::volScalarField& Foam::solvers::conductingFluid::getTemperature()
{   
    return thermo_.T();
}

// ************************************************************************* //
