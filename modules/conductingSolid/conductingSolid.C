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

#include "conductingSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(conductingSolid, 0);
    addToRunTimeSelectionTable(solver, conductingSolid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::conductingSolid::conductingSolid
(
    fvMesh& mesh
)
:
    solid(mesh),

    electro_
    (
        electromagneticModel::New(mesh)
    ),

    electro(electro_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::conductingSolid::~conductingSolid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


/*
void Foam::solvers::conductingSolid::setJJsigma(volScalarField& JJsigma)
{
    JJsigma_=JJsigma;
}
Foam::volScalarField& Foam::solvers::conductingSolid::getTemperature()//volScalarField& T_external)
{   
    return thermo_.T();
}
*/

void Foam::solvers::conductingSolid::postSolve()
{
    solid::postSolve();

    if (electro.correctElectromagnetics())
    {
        electromagneticPredictor();
    }
}

void Foam::solvers::conductingSolid::electromagneticPredictor()
{
    bool imaginary = electro.isComplex();
    //Get references for modification
    //Lorentz force term
    electro_->JxB =
    0.5*(
        (electro.J() ^ electro.B() )
        +(electro.J(imaginary) ^ electro.B(imaginary) )
    );
    //Joule heating
    //multiply by inverse of sigma to avoid division by zero
    electro_->JJsigma =
    0.5*(
        (electro.J() & electro.J())
        +(electro.J(imaginary) & electro.J(imaginary))
    )*electro.sigmaInv();

    electro_->setCorrected();
}

void Foam::solvers::conductingSolid::setCorrectElectromagnetics()
{
    electro_->setCorrectElectromagnetics();
}

// ************************************************************************* //
