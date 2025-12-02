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

#include "incompressibleConductingVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "findRefCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleConductingVoF, 0);
    addToRunTimeSelectionTable(solver, incompressibleConductingVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingVoF::incompressibleConductingVoF(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture
    (
        refCast<incompressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),
    electroBase(mesh),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),

    alpha1_old_
    (
        IOobject
        (
            "alpha1_old",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1
    ),

    alpha1_old(alpha1_old_),

    U_old_
    (
        IOobject
        (
            "U_old",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_
    ),

    U_old(U_old_),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        mixture
    )
{
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
    // Read the controls
    readControls();

    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingVoF::~incompressibleConductingVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getAlpha1()
{
    return alpha1;
}

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getAlpha2()
{
    return alpha2;
}

const Foam::word& Foam::solvers::incompressibleConductingVoF::getPhase1Name()
{
    return mixture.phase1Name();
}

const Foam::word& Foam::solvers::incompressibleConductingVoF::getPhase2Name()
{
    return mixture.phase2Name();
}

void Foam::solvers::incompressibleConductingVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    // Calculate the mass-flux from the accumulated alphaPhi1
    rhoPhi = (alphaPhi1*(rho1 - rho2) + phi*rho2);

    if (pimple.predictTransport())
    {
        momentumTransport.predict();
    }
}


void Foam::solvers::incompressibleConductingVoF::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::incompressibleConductingVoF::thermophysicalPredictor()
{}

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getAlpha1Old() const
{
    return alpha1_old_;
}

void Foam::solvers::incompressibleConductingVoF::storeAlpha1()
{
    alpha1_old_ = alpha1;
}

void Foam::solvers::incompressibleConductingVoF::storeU()
{
    U_old_ = U_;
}

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getNu()
{
    return mixture.nu();
}


void Foam::solvers::incompressibleConductingVoF::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport.correct();
    }
    if (electro.correctElectromagnetics())
    {
        // Set as corrected
        electro_.setCorrected();

        //Correct current density
        //electro_.correct();
    }
}

void Foam::solvers::incompressibleConductingVoF::solveElectromagnetics()
{
    // Doesn't solve for new currents, but at least makes sure that currents stay in the conducting phase.
    fixCurrentsWithMixtureConductivity();
    //Update JxB and Joule heating terms without solving
    electromagneticPredictor();
    //Solve potential equation
    //electro_.solve();
}


void Foam::solvers::incompressibleConductingVoF::fixCurrentsWithMixtureConductivity()
{
    volScalarField sigmaNormalized = electro.sigma()/
    dimensionedScalar
    (
        "sigmaMag",
        pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
        gMax(electro.sigma()())
    );
    volVectorField new_Jre = electro.J()*sigmaNormalized;
    volVectorField& Jre = getJ();
    forAll(new_Jre, cellI)
    {
        Jre[cellI] = new_Jre[cellI];
    }
    Jre.correctBoundaryConditions();
    if (electro.isComplex())
    {
        volVectorField new_Jim = electro.J(true)*sigmaNormalized;
        volVectorField& Jim = getJ(true);
        forAll(new_Jim, cellI)
        {
            Jim[cellI] = new_Jim[cellI];
        }
        Jim.correctBoundaryConditions();
    }
}

// ************************************************************************* //
