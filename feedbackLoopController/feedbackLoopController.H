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
    // Control values
    Pair<scalar> target_value_;
    word controlType_;
    bool first_iteration_ = true;
    Pair<int> sign_ = Pair<int>(1,1);
    Pair<int> reset_counter_ = Pair<int>(0,0);
    // Auxiliary values
    Pair<scalar> previous_error_ = Pair<scalar>(0,0);
    Pair<scalar> control_value_;
    Pair<scalar> integral_ = Pair<scalar>(0,0);
    Pair<scalar> oldTime_ = Pair<scalar>(0,0);
    bool needs_update_ = true;
    // Coefficient values
    Pair<scalar> proportional_coeff_ = Pair<scalar>(1,1);
    Pair<scalar> differential_coeff_ = Pair<scalar>(0,0);
    Pair<scalar> integral_coeff_ = Pair<scalar>(0,0);
    Pair<scalar> max_error_ = Pair<scalar>(0.05,0.05);//5%
    Pair<scalar> reference_value_ = Pair<scalar>(1,1);
    Pair<scalar> max_value_ = Pair<scalar>(great,great);
    Pair<scalar> min_value_ = Pair<scalar>(-great,-great);

    /*void updateControlValues_(
        const Pair<scalar>& target_value,
        const word& controlType);

    void updateAuxiliaryValues_(
        const Pair<scalar>& previous_error,
        const Pair<scalar>& control_value,
        const Pair<scalar>& integral,
        const scalar& oldTime);*/

public:

    // Public Data


    //- Runtime type information
    TypeName("feedbackLoopController");


    // Constructors

        //- Construct from region mesh
        feedbackLoopController();
        feedbackLoopController(
            const Pair<scalar>& target_value,
            const word& controlType,
            Pair<scalar> control_value = Pair<scalar>(0,0));

        //- Disallow default bitwise copy construction
        feedbackLoopController(const feedbackLoopController&) = delete;


    //- Destructor
    virtual ~feedbackLoopController();


    // Member Functions

        //- 
        Pair<scalar> calculateCorrection(Pair<scalar> present_value, scalar newTime);
        word getControlType();
        bool needsUpdate(scalar error, label id);
        bool needsUpdate(Pair<scalar> error);
        bool needsUpdate();
        void setReference(Pair<scalar> reference_value);
        void setMaxError(Pair<scalar> max_error);
        void setMinMaxValue(Pair<scalar> min_value,Pair<scalar> max_value);
        Pair<scalar> getControlValues();
        //void setControlType(word controlType);

        void updateCoefficients(
            const Pair<scalar>& proportional_coeff,
            const Pair<scalar>& differential_coeff,
            const Pair<scalar>& integral_coeff);


    // Member Operators

        //- Copy control values
        feedbackLoopController& operator=(const feedbackLoopController& newController)
        {
            if (this == &newController)  // Self-assignment check
            return *this;

            //TODO: check if usefult to copy other values
            this->target_value_ = newController.target_value_;
            this->controlType_ = newController.controlType_;
            this->control_value_ = newController.control_value_;
            return *this;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
