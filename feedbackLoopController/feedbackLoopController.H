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
    Foam::feedbackLoopController

Description
    Module for handling region/material electromagnetic properties.

SourceFiles
    feedbackLoopController.C

\*---------------------------------------------------------------------------*/

#ifndef feedbackLoopController_H
#define feedbackLoopController_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class feedbackLoopController Declaration
\*---------------------------------------------------------------------------*/

//Base class for electromagnetic modules
class feedbackLoopController
{

protected:
    Pair<scalar> target_value_;
    word controlType_;

    Pair<scalar> proportional_coeff_;
    Pair<scalar> differential_coeff_;
    Pair<scalar> integral_coeff_;

    Pair<scalar> previous_error_;
    Pair<scalar> control_value_;
    Pair<scalar> integral_;
    scalar oldTime_ = 0;

public:

    // Public Data


    //- Runtime type information
    TypeName("feedbackLoopController");


    // Constructors

        //- Construct from region mesh
        feedbackLoopController(
            const Pair<scalar>& target_value,
            const word& controlType);

        //- Disallow default bitwise copy construction
        feedbackLoopController(const feedbackLoopController&) = delete;


    //- Destructor
    virtual ~feedbackLoopController();


    // Member Functions

        //- 
        Pair<scalar> calculateCorrection(Pair<scalar> present_value, scalar newTime);
        word getControlType();
        //void setControlType(word controlType);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const feedbackLoopController&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
