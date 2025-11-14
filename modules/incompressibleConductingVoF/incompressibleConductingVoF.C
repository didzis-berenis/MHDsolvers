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
    twoPhaseConductingVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture
    (
        refCast<incompressibleTwoPhaseVoFMixture>(twoPhaseConductingVoFSolver::mixture)
    ),

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
    Info << "Getting alpha1 field" << endl;
    return mixture.alpha1();
}

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getAlpha2()
{
    Info << "Getting alpha2 field" << endl;
    return mixture.alpha2();
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
    twoPhaseConductingVoFSolver::prePredictor();

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


void Foam::solvers::incompressibleConductingVoF::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport.correct();
    }
    if (electro.correctElectromagnetics())
    {
        //Correct current density
        //electro_.correct();

        // Set as corrected
        electro_.setCorrected();
    }
}

void Foam::solvers::incompressibleConductingVoF::solveElectromagnetics()
{
    volScalarField sigmaNormalized = electro.sigma()/dimensionedScalar
    (
        "sigmaMag",
        pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
        gMax((electro.sigma())())
    );
    volVectorField new_Jre = electro.J()*sigmaNormalized;
    volVectorField& Jre = getJ();
    forAll(new_Jre, cellI)
    {
        Jre[cellI] = new_Jre[cellI];
    }
    Jre.correctBoundaryConditions();
    //Jre.write();
    if (electro.isComplex())
    {
        volVectorField new_Jim = electro.J(true)*sigmaNormalized;
        volVectorField& Jim = getJ(true);
        forAll(new_Jim, cellI)
        {
            Jim[cellI] = new_Jim[cellI];
        }
        Jim.correctBoundaryConditions();
        //Jim.write();
    }
    electromagneticPredictor();
    //Solve potential equation
    //electro_.solve();

    // Doesn't solve for new currents, but at least makes sure that currents stay in the conducting phase.
    /*volScalarField sigmaNormalized = electro.sigma()/gMax((electro.sigma())());
    volVectorField Jre = electro.J()*sigmaNormalized;
    electro_.setJ(Jre);
    if (electro.isHarmonic())
    {
        volVectorField Jim = electro.J(true)*sigmaNormalized;
        electro_.setJ(Jim,true);
    }
    electromagneticPredictor();*/
}


/*void Foam::solvers::incompressibleConductingVoF::postSolve()
{
    twoPhaseConductingSolver::postSolve();
    Pout << "Solving electromagnetics" << endl;
    volScalarField sigmaNormalized = electro.sigma()/dimensionedScalar
    (
        "sigmaMag",
        pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
        gMax((electro.sigma())())
    );
    volVectorField new_Jre = electro.J()*sigmaNormalized;
    volVectorField& Jre = getJ();
    forAll(new_Jre, cellI)
    {
        Jre[cellI] = new_Jre[cellI];
    }
    Jre.correctBoundaryConditions();
    //Jre.write();
    if (electro.isComplex())
    {
        volVectorField new_Jim = electro.J(true)*sigmaNormalized;
        volVectorField& Jim = getJ(true);
        forAll(new_Jim, cellI)
        {
            Jim[cellI] = new_Jim[cellI];
        }
        Jim.correctBoundaryConditions();
        //Jim.write();
    }
    electromagneticPredictor();
}*/

// ************************************************************************* //
