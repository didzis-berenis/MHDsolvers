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
void Foam::solvers::incompressibleConductingFluid::deltaJ(volVectorField& deltaU, bool imaginary)
{
    //Interpolating cross product u x B over mesh faces
    surfaceScalarField psiUB = fvc::interpolate(deltaU ^ electro.B(imaginary)) & mesh.Sf();
    //Get reference for modification
    volScalarField& PotE = electroPtr->PotE(imaginary);
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
        volVectorField deltaU = U_-U_old_;
        bool imaginary = electro.isComplex();
        //Get J difference by incorporating deltaU x B term
        volVectorField deltaJre = electro.deltaJ(deltaU);
        volVectorField deltaJim = deltaJre;
        if (imaginary)
        {
            deltaJim = electro.deltaJ(deltaU, imaginary);
        }

        //Get references for modification
        //Lorentz force term
        electroPtr->JxB =
        0.5*(
            ((electro.J()+deltaJre) ^ electro.B() )
            +((electro.J(imaginary)+deltaJim) ^ electro.B(imaginary) )
        );
        //Joule heating
        //multiply by inverse of sigma to avoid division by zero
        electroPtr->JJsigma =
        0.5*(
            ((electro.J()+deltaJre) & (electro.J()+deltaJre))
            +((electro.J(imaginary)+deltaJim) & (electro.J(imaginary)+deltaJim))
        )*electro.sigmaInv();
}


// ************************************************************************* //
