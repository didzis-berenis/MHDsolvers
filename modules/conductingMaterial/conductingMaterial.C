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

#include "conductingMaterial.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(conductingMaterial, 0);
    addToRunTimeSelectionTable(solver, conductingMaterial, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::conductingMaterial::conductingMaterial
(
    fvMesh& mesh
)
:
    solver(mesh),
    electroBase(mesh)
{
    // Non-source conductingMaterial regions are considered passive.
    // Initialize only source regions.
    //if (electro.isSource())
    //{
        label PotERefCell = 0;
        scalar PotERefValue = 0.0;
        setRefCell
        (
            electro.PotE(),
            pimple.dict(),
            PotERefCell,
            PotERefValue
        );
        mesh.schemes().setFluxRequired(electro.PotE().name());

        if (electro.isComplex())
        {
            setRefCell
            (
                electro.PotE(true),
                pimple.dict(),
                PotERefCell,
                PotERefValue
            );
            mesh.schemes().setFluxRequired(electro.PotE(true).name());
        }
    //}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::conductingMaterial::~conductingMaterial()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalar Foam::solvers::conductingMaterial::maxDeltaT() const
{
    scalar deltaT = fvModels().maxDeltaT();
    return deltaT;
}


void Foam::solvers::conductingMaterial::preSolve()
{}


void Foam::solvers::conductingMaterial::moveMesh()
{}


void Foam::solvers::conductingMaterial::prePredictor()
{}


void Foam::solvers::conductingMaterial::momentumPredictor()
{}


void Foam::solvers::conductingMaterial::thermophysicalPredictor()
{}


void Foam::solvers::conductingMaterial::pressureCorrector()
{}


void Foam::solvers::conductingMaterial::postCorrector()
{
    // Solve only for source regions
    if (electro.correctElectromagnetics())
    {
        //Calculate current density
        /*electro_.findJ();
        bool imaginary = electro.isComplex();
        if (imaginary)
        {
            electro_.findJ(imaginary);
        }
        if (electro.getRegionRole() == "coil")
        {
            electro_.predict();
        }
        else
        {
            const volVectorField Jre = electro.J();
            volVectorField& JreRef = electro_.Jref();
            JreRef = Jre;
            if (imaginary)
            {
                const volVectorField Jim = electro.J(imaginary);
                volVectorField& JimRef = electro_.Jref(imaginary);
                JimRef = Jim;
            }
        }*/
        electro_.setCorrected();
    }
}

void Foam::solvers::conductingMaterial::postSolve()
{}

void Foam::solvers::conductingMaterial::solveElectromagnetics()
{
    // Solve only for source regions
    //if (electro.isSource())
    //{
        //electro_.solve();
    //}
}

// ************************************************************************* //
