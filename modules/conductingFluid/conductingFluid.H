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
    Foam::solvers::fluid

Description
    Solver module for steady or transient turbulent flow of compressible fluids
    with heat-transfer for HVAC and similar applications, with optional
    mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian particles,
    radiation, surface film etc. and constraining or limiting the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    conductingFluid.C

\*---------------------------------------------------------------------------*/

#ifndef conductingFluid_H
#define conductingFluid_H

#include "isothermalFluid.H"
#include "electromagneticModel.H"
#include "fluidThermoThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class conductingFluid Declaration
\*---------------------------------------------------------------------------*/

class conductingFluid
:
    public isothermalFluid
{

protected:

    // Thermophysical transport

        autoPtr<fluidThermoThermophysicalTransportModel>
            thermophysicalTransport;

    //- Pointer to the electromagnetic model
        autoPtr<electromagneticModel> electroPtr;

    //- Store old velocity field
        volVectorField U_old_;

    //- Current correction
        //real
        //tmp<volVectorField> JUBre_;
        //imaginary
        //tmp<volVectorField> JUBim_;

        // Get velocity from solver
        //virtual tmp<volVectorField>& JUB(bool imaginary=false);

public:
        

    //- Runtime type information
    TypeName("conductingFluid");


    // Constructors

        //- Construct from region mesh
        conductingFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        conductingFluid(const conductingFluid&) = delete;


    //- Destructor
    virtual ~conductingFluid();


    // Member Functions

        // Assign volume force to the solver
        //virtual void setJJsigma(volScalarField&);

        // Assign volume heating to the solver
        //virtual void setJxB(volVectorField&);

        // Get pressure from solver
        //virtual volScalarField& getPressure();

        // Get temperature from solver
        //virtual volScalarField& getTemperature();

        //- Reference to the fluid fluid electromagnetic properties
        const electromagneticModel& electro;

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void momentumPredictor();
        virtual void thermophysicalPredictor();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Solve electromagnetic modelling
        virtual void postSolve();
        virtual void solveElectromagnetics();
        //virtual void correctJ(bool imaginary=false);
        //virtual void electromagneticPredictor();

        //- Set correctElectromagnetics switch to true
        virtual void setCorrectElectromagnetics();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const conductingFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
