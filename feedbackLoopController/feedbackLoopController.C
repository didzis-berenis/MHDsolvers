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

Foam::feedbackLoopController::feedbackLoopController()
:
    target_value_(0,0),
    controlType_("")
{}

Foam::feedbackLoopController::feedbackLoopController
(
    const Pair<scalar>& target_value,
    const word& controlType,
    Pair<scalar> control_value
)
:
    target_value_(target_value),
    controlType_(controlType),
    control_value_(control_value)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::feedbackLoopController::~feedbackLoopController()
{}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

/*void Foam::feedbackLoopController::updateControlValues_(
    const Pair<scalar>& target_value,
    const word& controlType)
{}

void Foam::feedbackLoopController::updateAuxiliaryValues_(
    const Pair<scalar>& previous_error,
    const Pair<scalar>& control_value,
    const Pair<scalar>& integral,
    const scalar& oldTime)
{}*/

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::Pair<Foam::scalar> Foam::feedbackLoopController::calculateCorrection(Pair<scalar> present_value, scalar newTime)
{
    Pair<scalar> integral = integral_;
    Pair<scalar> control_value = control_value_;
    Pair<scalar> previous_error = previous_error_;
    forAll(present_value, i)
    {
        scalar error = target_value_[i]-present_value[i];
        scalar proportional_correction = proportional_coeff_[i]*error;
        scalar differential_correction = 0;
        scalar integral_correction = 0;

        scalar deltaT = newTime - oldTime_;
        if (deltaT > 0)
        {
            differential_correction = differential_coeff_[i]*(error-previous_error[i])/deltaT;
            integral[i] += error*deltaT;
            integral_correction = integral_coeff_[i]*integral_[i];
        }
        previous_error[i] = error;
        control_value[i] +=  proportional_correction + differential_correction + integral_correction;
    }
    if (needsUpdate(previous_error))
    {
        integral_ = integral;
        previous_error_ = previous_error;
        control_value_ = control_value;
        oldTime_ = newTime;
    }
    return control_value_;
}

Foam::word Foam::feedbackLoopController::getControlType()
{
    return controlType_;
}
void Foam::feedbackLoopController::setReference(Pair<scalar> reference_value)
{
    forAll(reference_value, i)
    {
        if (reference_value[i] != 0)
            reference_value_[i] = reference_value[i];
    }
}
void Foam::feedbackLoopController::setMaxError(Pair<scalar> max_error)
{
    max_error_ = max_error;
}

bool Foam::feedbackLoopController::needsUpdate(Pair<scalar> error)
{
    needs_update_ = (std::abs(error[0]) > std::abs(max_error_[0]*reference_value_[0]) ||
        std::abs(error[1]) > std::abs(max_error_[1]*reference_value_[1]));
    return needsUpdate();
}

bool Foam::feedbackLoopController::needsUpdate()
{
    return needs_update_;
}

Foam::Pair<Foam::scalar> Foam::feedbackLoopController::getControlValues()
{
    return control_value_;
}

void Foam::feedbackLoopController::updateCoefficients(
    const Pair<scalar>& proportional_coeff,
    const Pair<scalar>& differential_coeff,
    const Pair<scalar>& integral_coeff)
{
    proportional_coeff_ = proportional_coeff;
    differential_coeff_ = differential_coeff;
    integral_coeff_ = integral_coeff;
}

// ************************************************************************* //
