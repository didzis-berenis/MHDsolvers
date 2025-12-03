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
#include "feedbackLoopController.H"
#include <map>
#include <vector>
#include <fstream>

static const double MAX_CONTROL_ERROR = 0.05;//relative error of 5%

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

        // Map from terminal name to wire regions
        std::map<word, std::vector<word>> terminalToRegions_;

        //- Feedback loop controllers
        std::map<word, bool> voltageNeedsUpdate_;
        std::map<word, feedbackLoopController> feedbackControllers_;
        std::map<word, scalar> inducedVoltageValue_;
        std::map<word, scalar> inducedVoltagePhase_;
        struct coilParameters_
        {
            word terminal_name;
            word region_name;
            scalar target_value;
            scalar target_phase;
            scalar value_multiplier;
            Foam::vector winding_direction;
        };
        scalar target_value_error_;
        scalar target_phase_error_;
        scalar induced_voltage_relaxation_factor_;
        std::map<word, coilParameters_> voltageControls_;
        std::map<word, vectorField> regionOldB_;
        std::map<word, std::vector<word>> terminalCoreRegions_;

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
        long long int cleanupCounter_ = 0;
        int restartInterval_;
        int integrationSteps_ = 0;
        int integrationCounter_ = 0;
        bool controllerCheckTime_ = false;
        std::map<word, scalar> integralCurrent_;
        scalar tStepElectro_ = 0;

        // Electromagnetic excitation frequency
        const scalar frequency_;

        const word fluidSolverName_ = "conductingFluid";
        const word solidSolverName_ = "conductingSolid";
        const word electricSolverName_ = "conductingMaterial";
        const word magneticSolverName_ = "magneticMaterial";
        //Perhaps better to construct set, if more solvers are implemented.
        //std::set<word> conductingFluidNames_;
        //std::set<word> conductingSolidNames_;

        //- Number of potential corrections
        int nPotEcorr_ = 3;//Optimal value for minimizing error and avoiding too many PIMPLE iterations
        int PotEcorr_ = 0;
        //- Get pointer electroBase class from region
        electroBase* getElectroBasePtr_(const word regionName);

        //- Get pointers to region solvers
        solvers::conductingFluid* getFluidPtr_(const word regionName);
        solvers::conductingSolid* getSolidPtr_(const word regionName);

        // Map global cell values to region cell values
        void scalarGlobalToField_(volScalarField& global,volScalarField& region,const word& regionName);
        void vectorGlobalToField_(volVectorField& global,volVectorField& region,const word& regionName);

        //- Global mesh pointer
        autoPtr<fvMesh> globalMesh_;
        //- Construct new global mesh from region meshes
        autoPtr<fvMesh> globalMeshNew_();

        //- Initialize electric potential boundary conditions
        void evaluatePotEBfs_(const word regionName, bool imaginary=false);
        void updatePotErefGrad_(const word regionName, scalar newGrad, bool imaginary=false);
        void updatePotErefValue_(const word terminalName, scalar newValue, bool imaginary=false);
        scalar getPotErefValue_(const word terminalName, bool imaginary=false);
        scalar calculateBoundaryCurrent_(const word terminalName, bool imaginary=false);
        int initial_iter_ = 0;
        int wait_iter_ = 0;

        //- Initialize current correction boundary conditions
        void evaluateDeltaJBfs_(const word regionName, bool imaginary=false);
        void evaluateJBfs_(const word regionName, bool imaginary=false);

        bool hasVoltageControl_=false;
        void checkIfAnyVoltageControl_();
        //void checkIfAnyElectricSources_()

        // Performs feedback control initialization
        void setUpFeedbackControllers_();
        int surfaceHitCount_(const point otherRegionPoint,const Foam::vector perpendicularVector,const fvPatch&);

        // Calculates phase shift for transient simulation
        scalar getTransientPhaseShift_(scalar argument, scalar step_shift);
        // Calculates unit vector field from current density
        volVectorField getJdirection_(word regionName, bool imaginary);
        //Check if imaginary Jref is the same direction as real part
        bool isJsameDirection(word regionName);
        // Calculates integral current density for a region
        scalar getCurrentSum_(word regionName,bool imaginary=false);//volVectorField JGlobal,
        scalar getInductionSum_(word regionName, coilParameters_ terminalControls, bool imaginary=false);
        // Writes provided value to a file
        void writeControlValue_(const word fileName, const scalar outputValue);
        // Reads value from a file
        scalar readControlValue_(const word fileName);
        bool fileExists_(const word fileName);

public:
    //constants
        const dimensionedScalar smallU = dimensionedScalar("smallU",dimVelocity,1e-6);
        static constexpr double PI = 3.141592653589793;
        static constexpr scalar mu_0 = 4*PI*(std::pow(10,-7));
        //- Magnetic field update controls
        const int waitInterval;

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

        //- Return region fvMesh
        fvMesh& mesh(const word regionName);

        scalar getArea(const word regionName, const word patchName);

        Foam::vector getCenter(const word regionName, const word patchName);

        Foam::vector getCenter(const word regionName);

        Foam::vector getNormal(const word regionName, const word patchName);

        void updateFeedbackControl();
        bool controllersNeedUpdate();

        // Return global fvMesh
        const fvMesh& globalMesh();

        //- Map for region to global cell id
        std::map<std::pair<word, label>, label> regionToGlobalCellId;

        //- Return region electromagneticModel
        const Foam::electromagneticModel& getElectro(const word regionName);

        //- Map region cell values to global cell values for scalar field
        void scalarFieldToGlobal(volScalarField& global,const volScalarField& region,const word& regionName);

        //- Map region cell values to global cell values for a vector field
        void vectorFieldToGlobal(volVectorField& global,const volVectorField& region,const word& regionName);

        //- Assign current density to solvers
        void setJ(volVectorField& field, bool imaginary=false);

        // Check if any region is electric source
        bool hasElectricSources();

        // Get formatted current density from electric source regions
        //volVectorField getSourceJ(word regionName, bool imaginary);

        //- Assign current density to specified region
        void setJToRegion(volVectorField& field, const word& regionName, bool imaginary=false);
        void setJRefToRegion(volVectorField& field, const word& regionName, bool imaginary=false);
        void setJRef(const word& regionName);

        //- Assign magnetic field to solvers
        void setB(volVectorField& field, bool imaginary=false);

        //- save old velocity
        void storeU(const word regionName);

        //- Check if solver type is fluidSolverName_
        bool isFluid(const word regionName);

        //- Check if solver type is solidSolverName_
        bool isSolid(const word regionName);

        //- Check if solver type is electricSolverName_
        bool isElectric(const word regionName);

        //- Check if solver type is magneticSolverName_
        bool isMagnetic(const word regionName);

        // Check if region is electric source
        bool needsReference(const word regionName);
        bool referenceNotSet(const word regionName);
        
        bool hasAnyRole(const word regionRole);
        
        //- Check if any region has velocity field
        bool checkIfAnyVelocityRegions();

        //- Check if electromagnetic solver time type is harmonic
        bool isElectroHarmonic();

        //- Check if magnetic field needs updating
        bool updateMagneticField();

        //- solve electromagnetics for region
        void solveElectromagnetics(const word regionName);

        //- Update electromagnetic system without solving
        void electromagneticPredictor(const word regionName);

        //- Set number of potential correctors
        void setPotentialCorrectors(const dictionary& dict);

        //- Check if potential corrector is needed
        bool correctElectroPotential();

        //- Initialize temperature gradient
        void calcTemperatureGradient(const word regionName);

        //- Check if cleanup of time directories is needed
        bool needsCleanup();

        //- Increment counter, which controls cleanup
        void countToCleanup();

        //- Return fluid solver
        const Foam::solvers::conductingFluid& getFluid(const word regionName);

        //- Return solid solver
        const Foam::solvers::conductingSolid& getSolid(const word regionName);

        //- Get region-solver name pairs
        List<Pair<word>> getNames();

        // Get solver id from region name
        int getRegionId(const word regionName);


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
