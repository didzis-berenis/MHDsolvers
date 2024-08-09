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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::conductingFluid::electromagneticPredictor()
{
    //Interpolating cross product u x B over mesh faces
    surfaceScalarField psiUB = fvc::interpolate((U_-U_old_) ^ electro.B()) & mesh.Sf();
    volScalarField& PotE = electro_->PotE();
    volScalarField sigma = electro.sigma();
    //Poisson equation for electric potential
    fvScalarMatrix PotEEqn
    (
    fvm::laplacian(sigma,PotE)
    ==
    sigma*fvc::div(psiUB)
    );
    //Reference potential
    label PotERefCell = 0;
    scalar PotERefValue = 0.0;
    PotEEqn.setReference(PotERefCell, PotERefValue);
    //Solving Poisson equation
    PotEEqn.solve();

    //Computation of current density at cell faces
    surfaceScalarField En = -(fvc::snGrad(PotE) * mesh.magSf()) + psiUB;
    //Current density at face center
    surfaceVectorField Env = En * mesh.Cf();
    
    //Get boundary conditions from J
    volVectorField JUB = electro.J();
    //Interpolation of current density at cell center
    JUB = sigma*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh.C()) );
    //Update current density distribution and boundary condition
    JUB.correctBoundaryConditions();

    //Lorentz force term
    electro_->JxB =  ((electro.J()+JUB) ^ electro.B() );
    //Joule heating
	electro_->JJsigma = ((electro.J()+JUB) & (electro.J()+JUB))*electro.sigmaInv();//multiply by inverse of sigma

    U_old_ = U_;
}


// ************************************************************************* //
