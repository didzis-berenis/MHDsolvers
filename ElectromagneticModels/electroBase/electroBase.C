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

#include "electroBase.H"
//#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroBase, 0);
    //defineRunTimeSelectionTable(electroBase, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroBase::electroBase
(
    fvMesh& mesh
)
:
    electroPtr_(electromagneticModel::New(mesh)),
    electro_(electroPtr_()),
    electro(electroPtr_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroBase::~electroBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::electroBase::solveElectromagnetics()
{
    if (electro.correctElectromagnetics())
    {
        //Correct current density
        electroPtr_->correct();
    }
}

void Foam::electroBase::electromagneticPredictor()
{
    electroPtr_->predict();
}

void Foam::electroBase::setCorrectElectromagnetics()
{
    electroPtr_->setCorrectElectromagnetics();
}

Foam::volVectorField& Foam::electroBase::getJ(bool imaginary)
{
    return electroPtr_->J(imaginary);
}

Foam::volVectorField& Foam::electroBase::getB(bool imaginary)
{
    return electroPtr_->B(imaginary);
}

// ************************************************************************* //
