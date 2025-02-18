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

#include "feedbackLoopController.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(feedbackLoopController, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::feedbackLoopController::feedbackLoopController
(
    const Pair<scalar>& target_value,
    const word& controlType
)
:
    target_value_(target_value),
    controlType_(controlType),
    proportional_coeff_(1,1),
    differential_coeff_(0,0),
    integral_coeff_(0,0),
    previous_error_(0,0),
    control_value_(0,0),
    integral_(0,0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::feedbackLoopController::~feedbackLoopController()
{}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::Pair<Foam::scalar> Foam::feedbackLoopController::calculateCorrection(Pair<scalar> present_value, scalar newTime)
{
    forAll(present_value, i)
    {
        scalar error = target_value_[i]-present_value[i];
        scalar proportional_correction = proportional_coeff_[i]*error;
        scalar differential_correction = 0;
        scalar integral_correction = 0;

        scalar deltaT = newTime - oldTime_;
        if (deltaT > 0)
        {
            differential_correction = differential_coeff_[i]*(error-previous_error_[i])/deltaT;
            integral_[i] += error*deltaT;
            integral_correction = integral_coeff_[i]*integral_[i];
        }
        previous_error_[i] = error;
        control_value_[i] =  proportional_correction + differential_correction + integral_correction;
    }
    oldTime_ = newTime;
    return control_value_;
}

Foam::word Foam::feedbackLoopController::getControlType()
{
    return controlType_;
}

// ************************************************************************* //
