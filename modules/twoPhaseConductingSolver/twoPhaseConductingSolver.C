/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "twoPhaseConductingSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(twoPhaseConductingSolver, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::twoPhaseConductingSolver::twoPhaseConductingSolver
(
    fvMesh& mesh,
    autoPtr<twoPhaseVoFMixture> mixturePtr
)
:
    conductingVoFSolver(mesh, autoPtr<VoFMixture>(mixturePtr.ptr())),

    mixture(refCast<twoPhaseVoFMixture>(conductingVoFSolver::mixture_)),

    alpha1(mixture.alpha1()),
    alpha2(mixture.alpha2()),

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

    alphaRestart
    (
        typeIOobject<surfaceScalarField>
        (
            IOobject::groupName("alphaPhi", alpha1.group()),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ).headerOk()
    ),

    alphaPhi1
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", alpha1.group()),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*fvc::interpolate(alpha1)
    )
{
    mesh.schemes().setFluxRequired(alpha1.name());

    if (alphaRestart)
    {
        Info << "Restarting alpha" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::twoPhaseConductingSolver::~twoPhaseConductingSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::twoPhaseConductingSolver::preSolve()
{
    conductingVoFSolver::preSolve();

    // Do not apply previous time-step mesh compression flux
    // if the mesh topology changed
    if (mesh().topoChanged())
    {
        talphaPhi1Corr0.clear();
    }
}


void Foam::solvers::twoPhaseConductingSolver::prePredictor()
{
    conductingVoFSolver::prePredictor();
    alphaPredictor();
    mixture.correct();
}

const Foam::volScalarField& Foam::solvers::twoPhaseConductingSolver::getAlpha1Old() const
{
    return alpha1_old_;
}

void Foam::solvers::twoPhaseConductingSolver::storeAlpha1()
{
    alpha1_old_ = alpha1;
}

void Foam::solvers::twoPhaseConductingSolver::storeU()
{
    U_old_ = U_;
}

void Foam::solvers::twoPhaseConductingSolver::momentumPredictor()
{
    conductingVoFSolver::momentumPredictor();
    storeAlpha1();
    storeU();
    // Might be unnecessary since electromagnetics is updated more often, due to surface change
    // Update deltaU for electromagnetic model
    volVectorField deltaU = U_ - U_old_;
    electro_.updateDeltaU(deltaU);
}


// ************************************************************************* //
