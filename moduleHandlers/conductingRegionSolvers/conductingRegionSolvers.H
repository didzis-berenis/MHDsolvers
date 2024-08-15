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
    Foam::conductingRegionSolvers

Description
    Class to hold the lists of region meshes and solvers

    Also provides loop and iteration functionality which automatically set the
    region Info prefix for each of the solvers returned by the '[]' operator or
    iterator.

Usage
    Given the \c conductingRegionSolvers named solvers:
    \verbatim
        // Create the region meshes and solvers
        conductingRegionSolvers solvers(runTime);
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
        forAllIter(conductingRegionSolvers, solvers, solver)
        {
            solver->momentumPredictor();
        }
    \endverbatim
    where the iterator increment sets the region \c Info prefix and at the end
    automatically resets the \c Info prefix to spaces.

SourceFiles
    conductingRegionSolvers.C

\*---------------------------------------------------------------------------*/

#ifndef conductingRegionSolvers_H
#define conductingRegionSolvers_H

#include "solver.H"
#include "conductingFluid.H"
#include "conductingSolid.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class conductingRegionSolvers Declaration
\*---------------------------------------------------------------------------*/

class conductingRegionSolvers
{
    // Private Member Data
        const Time& runTime_;

        //- Region meshes
        PtrList<fvMesh> regions_;

        //- Region solvers
        PtrList<solver> solvers_;

        //- Global space padding prefix
        string prefix0_;

        //- List of space padded region prefixes
        stringList prefixes_;

        //- Region name and solver name pairs
        List<Pair<word>> names_;

        std::map<word, label> regionIdx_;

        //- Region characteristic geometric size
        List<scalar> characteristicSizes_;
        
        //- Magnetic field update controls
        scalar maxRemDiff_ = 0;
        scalar maxRelDiff_ = 1;

        //- Write controls
        int writeMultiplier_ = 1;
        long long int writeCounter_ = 0;
        word writeControlDict_ = word::null;
        bool adjustableRunTime_ = false;

        const word fluidSolverName_ = "conductingFluid";
        const word solidSolverName_ = "conductingSolid";
        //Perhaps better to construct set, if more solvers are implemented.
        //std::set<word> conductingFluidNames_;
        //std::set<word> conductingSolidNames_;

        electroBase* getElectroBasePtr_(const word regionName);
        solvers::conductingFluid* getFluidPtr_(const word regionName);
        solvers::conductingSolid* getSolidPtr_(const word regionName);
        void scalarGlobalToField_(volScalarField& global,volScalarField& region,const word& regionName);
        void vectorGlobalToField_(volVectorField& global,volVectorField& region,const word& regionName);

        //- Region solvers
        autoPtr<fvMesh> globalMesh_;


public:
    //constants
        const dimensionedScalar smallU = dimensionedScalar("smallU",dimVelocity,1e-6);
        static constexpr scalar mu_0 = 4*3.14159*(std::pow(10,-7));

    // Constructors

        //- Construct from components
        conductingRegionSolvers(const Time& runTime);

        //- Disallow default bitwise copy construction
        conductingRegionSolvers(const conductingRegionSolvers&) = delete;


    //- Destructor
    ~conductingRegionSolvers();


    // Member Functions

        //- Return the number of region solvers
        inline label size() const;

        //- Set the Info prefix to space padding for global messages
        void setGlobalPrefix() const;

        //- Set the Info prefix to the space padded region name
        void setPrefix(const label i) const;

        //- Reset the Info prefix to null
        void resetPrefix() const;

        // Return region fvMesh
        fvMesh& mesh(const word regionName);

        const fvMesh& globalMesh();

        std::map<std::pair<word, label>, label> localToGlobalID;

        void scalarFieldToGlobal(volScalarField& global,const volScalarField& region,const word& regionName);
        void vectorFieldToGlobal(volVectorField& global,const volVectorField& region,const word& regionName);
        const Foam::electromagneticModel& getElectro(const word regionName);

        // Assign current density to solvers
        void setJ(volVectorField& field, bool imaginary=false);

        // Assign magnetic field to solvers
        void setB(volVectorField& field, bool imaginary=false);

        // save old velocity
        void storeU(const word regionName);

        // Check if solver type is fluidSolverName_
        bool isFluid(const word regionName);

        // Check if solver type is solidSolverName_
        bool isSolid(const word regionName);

        // Check if electromagnetic solver time type is harmonic
        bool isElectroHarmonic();

        // Check if magnetic field needs updating
        bool updateMagneticField();

        //- Set correctElectromagnetics switch to true for all conducting regions
        void setCorrectElectromagnetics();

        //- solve electromagnetics for region
        void solveElectromagnetics(const word regionName);

        //- Update electromagnetic system without solving
        void electromagneticPredictor(const word regionName);
        
        const Foam::solvers::conductingFluid& getFluid(const word regionName);
        
        const Foam::solvers::conductingSolid& getSolid(const word regionName);

        // Get names
        List<Pair<word>> getNames();


    // Iterator

        class iterator;
        friend class iterator;

        // Iterator class to loop over the solvers, setting the prefix for each
        class iterator
        {
            // Private Data

                //- Reference to the conductingRegionSolvers
                conductingRegionSolvers& conductingRegionSolvers_;

                //- Current solver index
                label index_;

        public:

            friend class conductingRegionSolvers;

            // Constructors

                //- Construct for the conductingRegionSolvers
                inline explicit iterator(conductingRegionSolvers&);

            // Member Operators

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline solver& operator*();
                inline solver& operator()();
                inline solver* operator->();

                inline iterator operator++();
                inline iterator operator++(int);
        };

        //- Return an iterator to begin traversing the solvers
        inline iterator begin();

        //- Return an iterator to end traversing the solvers
        inline iterator end();


    // Member Operators

        //- Cast to the list of solvers
        inline operator PtrList<solver>&();

        //- Set the region i prefix and return the corresponding solver
        solver& operator[](const label i);
        solver& operator()(const word regionName);
        //solver* operator()(const label i);

        //- Disallow default bitwise assignment
        void operator=(const conductingRegionSolvers&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "conductingRegionSolversI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
