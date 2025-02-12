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
        const Time& runTime_;

        //- Solver pointer
        solver* solverPtr_;

        //- Solver name
        word name_;

        //- Region characteristic geometric size
        scalar characteristicSize_;
        
        //- Magnetic field update controls
        scalar maxRemDiff_ = 0;
        scalar maxRelDiff_ = 1;

        //- Write controls
        int writeMultiplier_ = 1;
        long long int writeCounter_ = 0;
        word writeControlDict_ = word::null;
        bool adjustableRunTime_ = false;
        long long int cleanupCounter_ = 0;
        int restartInterval_;

        //- Names of derived solvers
        const word incompressibleFluidSolverName_ = "incompressibleConductingFluid";
        const word fluidSolverName_ = "conductingFluid";
        const word solidSolverName_ = "conductingSolid";

        //- Get pointer electroBase class from region
        electroBase* getElectroBasePtr_();

        //- Get pointers to derived solver class object
        solvers::incompressibleConductingFluid* getIncompressibleFluidPtr_();
        solvers::conductingFluid* getFluidPtr_();
        solvers::conductingSolid* getSolidPtr_();


public:
    //constants
        const dimensionedScalar smallU = dimensionedScalar("smallU",dimVelocity,1e-6);
        static constexpr double PI = 3.141592653589793;
        static constexpr scalar mu_0 = 4*PI*(std::pow(10,-7));
        //- Magnetic field update controls
        const int waitInterval;

    // Constructors

        //- Construct from components
        conductingRegionSolver(const Time& runTime, fvMesh& mesh);

        //- Disallow default bitwise copy construction
        conductingRegionSolver(const conductingRegionSolver&) = delete;


    //- Destructor
    ~conductingRegionSolver();


    // Member Functions
        const Foam::electromagneticModel& getElectro();

        // Assign current density to solvers
        void setJ(volVectorField& field, bool imaginary=false);

        // Assign magnetic field to solvers
        void setB(volVectorField& field, bool imaginary=false);

        // save old velocity
        void storeU();

        // Check if solver type is incompressibleFluidSolverName_
        bool isIncompressibleFluid();

        // Check if solver type is fluidSolverName_
        bool isFluid();

        // Check if solver type is solidSolverName_
        bool isSolid();

        // Check if electromagnetic solver time type is harmonic
        bool isElectroHarmonic();

        // Check if magnetic field needs updating
        bool updateMagneticField();

        //- solve electromagnetics for region
        void solveElectromagnetics();

        //- Update electromagnetic system without solving
        void electromagneticPredictor();

        //- Initialize temperature gradient
        void calcTemperatureGradient();

        //- Check if cleanup of time directories is needed
        bool needsCleanup();

        //- Increment counter, which controls cleanup
        void countToCleanup();

        // Read-only access

        // Get Velocity
        const volVectorField& getVelocity();

        // Get Temperature
        const volScalarField& getTemperature();

        //- Return incompressible fluid solver
        const Foam::solvers::incompressibleConductingFluid& getIncompressibleFluid();

        //- Return fluid solver
        const Foam::solvers::conductingFluid& getFluid();

        //- Return solid solver
        const Foam::solvers::conductingSolid& getSolid();


        // Getsolver name
        word getName();

        // Get reference to solver
        inline solver& operator()()
        {
            return *solverPtr_;//solver_();
        }

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
