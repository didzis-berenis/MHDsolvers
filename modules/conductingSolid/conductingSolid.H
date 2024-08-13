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

Class
    Foam::solvers::solid

Description
    Solver module for thermal transport in solid domains and regions for
    conjugate heat transfer, HVAC and similar applications, with optional mesh
    motion and mesh topology changes.

SourceFiles
    conductingSolid.C

\*---------------------------------------------------------------------------*/

#ifndef conductingSolid_H
#define conductingSolid_H

#include "solid.H"
#include "electromagneticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class conductingSolid Declaration
\*---------------------------------------------------------------------------*/

class conductingSolid
:
    public solid
{

protected:

    //- Pointer to the electromagnetic model
        autoPtr<electromagneticModel> electroPtr;

public:

    // Public Data


    //- Runtime type information
    TypeName("conductingSolid");


    // Constructors

        //- Construct from region mesh
        conductingSolid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        conductingSolid(const conductingSolid&) = delete;


    //- Destructor
    virtual ~conductingSolid();


    // Member Functions

        // Assign volume heating to the solver
        //virtual void setJJsigma(volScalarField&);

        // Get temperature from solver
        //virtual volScalarField& getTemperature();

        //- Reference to the fluid fluid electromagnetic properties
        const electromagneticModel& electro;
        virtual void postSolve();

        //- Set correctElectromagnetics switch to true
        virtual void setCorrectElectromagnetics();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const conductingSolid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
