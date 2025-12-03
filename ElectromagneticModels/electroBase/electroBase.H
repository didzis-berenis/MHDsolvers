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
    Foam::electroBase

Description
    Module for handling region/material electromagnetic properties.

SourceFiles
    electroBase.C

\*---------------------------------------------------------------------------*/

#ifndef electroBase_H
#define electroBase_H

#include "electromagneticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class electroBase Declaration
\*---------------------------------------------------------------------------*/

//Base class for electromagnetic modules
class electroBase
{
    //- Pointer to the electromagnetic model
        autoPtr<electromagneticModel> electroPtr_;

protected:

        //- Non-const reference to the fluid fluid electromagnetic properties
        electromagneticModel& electro_;

        bool current_reference_set_re_ = false;

        bool current_reference_set_im_ = false;

public:

    // Public Data


    //- Runtime type information
    TypeName("electroBase");


    // Constructors

        //- Construct from region mesh
        electroBase(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        electroBase(const electroBase&) = delete;


    //- Destructor
    virtual ~electroBase();


    // Member Functions

        //- Reference to the fluid fluid electromagnetic properties
        const electromagneticModel& electro;
        
        //- Solve electromagnetic system
        virtual void solveElectromagnetics();

        //- Update electromagnetic system without solving
        virtual void electromagneticPredictor();

        //- Electric conductivity [S/m]
        virtual volScalarField& getSigma();

        //- Current density [A/m^2]
        virtual volVectorField& getJ(bool imaginary=false);

        //- Current density reference value [A/m^2]
        virtual volVectorField& getJref(bool imaginary=false);
        //- Sets current_reference_set_ flag to true
        virtual void markCurrentReferenceAsSet(bool imaginary=false);
        //- Returns current_reference_set_ value
        virtual bool currentReferenceSet(bool imaginary=false);

        //- Magnetic field [T]
        virtual volVectorField& getB(bool imaginary=false);

        //virtual void evalSigmaBnd();
        scalar getI(const word terminalName, bool imaginary = false);
        //- Electric potential [V/m]
        virtual void initPotE(bool imaginary=false);
        virtual void updatePotErefGrad(scalar newGrad, bool imaginary);
        virtual void updatePotErefValue(const word terminalName, scalar newValue, bool imaginary);
        scalar getPotErefValue(const word terminalName, bool imaginary);
        bool hasBoundary(const word terminalName);
        bool isElectricPotentialBoundary(const word terminalName);
        word findGroundTerminal(const word terminalName);
        Foam::vector getCenter(const word terminalName);

        //- Current density correction [A/m^2]
        virtual void initDeltaJ(bool imaginary=false);
        virtual void initJ(bool imaginary=false);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const electroBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
