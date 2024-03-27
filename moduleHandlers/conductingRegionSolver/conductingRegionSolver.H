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
    Foam::conductingRegionSolver

Description
    Class to hold the lists of region meshes and solvers

    Also provides loop and iteration functionality which automatically set the
    region Info prefix for each of the solvers returned by the '[]' operator or
    iterator.

Usage
    Given the \c conductingRegionSolver named solvers:
    \verbatim
        // Create the region meshes and solvers
        conductingRegionSolver solvers(runTime);
    \endverbatim

    The list of solvers can be looped over:
    \verbatim
        forAll(solvers, i)
        {
            solvers[i].momentumPredictor();
        }
    \endverbatim
    where the '[]' operator sets the region Info prefix.  After the loop the
    region Info prefix remains set to the last region prefix and so for global
    messages, e.g. the global time-step the Info prefix must be specifically
    reset to spaces by calling the \c setGlobalPrefix() function.

    Alternatively the list of solvers can be iterated over:
    \verbatim
        forAllIter(conductingRegionSolver, solvers, solver)
        {
            solver->momentumPredictor();
        }
    \endverbatim
    where the iterator increment sets the region \c Info prefix and at the end
    automatically resets the \c Info prefix to spaces.

SourceFiles
    conductingRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef conductingRegionSolver_H
#define conductingRegionSolver_H

#include "solver.H"
#include "incompressibleConductingFluid.H"
#include "conductingFluid.H"
#include "conductingSolid.H"
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class conductingRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class conductingRegionSolver
{
    // Private Member Data

        //- Region solvers
        solver* solverPtr_;
        autoPtr<solver> solver_;

        //- Region name and solver name pairs
        word name_;

        const word incompressibleFluidSolverName_ = "incompressibleConductingFluid";
        const word fluidSolverName_ = "conductingFluid";
        const word solidSolverName_ = "conductingSolid";

        solvers::incompressibleConductingFluid* getIncompressibleFluidPtr_();

        solvers::conductingFluid* getFluidPtr_();

        solvers::conductingSolid* getSolidPtr_();


public:

    // Constructors

        //- Construct from components
        conductingRegionSolver(const Time& runTime, fvMesh& mesh);

        //- Disallow default bitwise copy construction
        conductingRegionSolver(const conductingRegionSolver&) = delete;


    //- Destructor
    ~conductingRegionSolver();


    // Member Functions

        // Return region fvMesh
        void setJxB(volVectorField& field);

        // Return region fvMesh
        void setJJsigma(volScalarField& field);

        // Return region fvMesh
        void getVelocity(volVectorField& field);

        // Return region fvMesh
        void getPressure(volScalarField& field);

        // Return region fvMesh
        void getTemperature(volScalarField& field);

        // Return region fvMesh
        bool isIncompressibleFluid();

        // Return region fvMesh
        bool isFluid();

        // Return region fvMesh
        bool isSolid();

        // Get names
        //word getName();

        inline solver& operator()();

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::solver& Foam::conductingRegionSolver::operator()()
{
    return solver_();
}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
