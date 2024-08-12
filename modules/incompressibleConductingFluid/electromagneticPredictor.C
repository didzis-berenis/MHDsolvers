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

#include "incompressibleConductingFluid.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::solvers::incompressibleConductingFluid::correctJ(bool imaginary)
{
    //Interpolating cross product u x B over mesh faces
    surfaceScalarField psiUB = fvc::interpolate((U_-U_old_) ^ electro.B(imaginary)) & mesh.Sf();
    //Get reference for modification
    volScalarField& PotE = electro_->PotE(imaginary);
    //Const access
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
    if (!JUB(imaginary).valid())
    {
        JUB(imaginary) = electro.J(imaginary);
    }
    volVectorField& JUBRef = JUB(imaginary).ref();
    //Interpolation of current density at cell center
    JUBRef = sigma*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh.C()) );
    //Update current density distribution and boundary condition
    JUBRef.correctBoundaryConditions();
}

void Foam::solvers::incompressibleConductingFluid::electromagneticPredictor()
{
        correctJ();
        bool imaginary = electro.isComplex();
        if (imaginary)
        {
            correctJ(imaginary);
        }

        //Get references for modification
        //Lorentz force term
        electro_->JxB =
        0.5*(
            ((electro.J()+JUB()) ^ electro.B() )
            +((electro.J(imaginary)+JUB(imaginary)) ^ electro.B(imaginary) )
        );
        //Joule heating
        //multiply by inverse of sigma to avoid division by zero
        electro_->JJsigma =
        0.5*(
            ((electro.J()+JUB()) & (electro.J()+JUB()))
            +((electro.J(imaginary)+JUB(imaginary)) & (electro.J(imaginary)+JUB(imaginary)))
        )*electro.sigmaInv();
        //Store old velocity for next update
        U_old_ = U_;
        electro_->setCorrected();
}


// ************************************************************************* //
