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
    Pair<scalar> oldTime = oldTime_;
    forAll(present_value, i)
    {
        scalar error = present_value[i]-target_value_[i]+offset_value_[i];
        // Switch sign to avoid positive feedback loop
        // This can happen if small error remains despite of incrementing control value
        if (!first_iteration_[i] && use_stabilizer_[i] && abs(error)>abs(previous_error[i]))
        {
            positive_feedback_counter_[i] ++;
        }
        if (positive_feedback_counter_[i]>positive_feedback_tolerance_[i])
        {
            offset_value_[i] = -error;
            Info << "Positive feedback loop detected!" << endl
            << "Setting present error as an offset value: offset_value_[" << i << "]: " 
            << -error << endl;
            // Reset to zero
            positive_feedback_counter_[i] = 0;
            first_iteration_[i] = true;
            error = present_value[i]-target_value_[i]+offset_value_[i];
        }
        scalar proportional_correction =
        proportional_coeff_[i]*error;
        scalar differential_correction = 0;
        scalar integral_correction = 0;

        scalar deltaT = newTime - oldTime[i];
        if (deltaT > 0)
        {
            differential_correction = differential_coeff_[i]*(error-previous_error[i])/deltaT;
        }
        integral[i] += error*deltaT;
        integral_correction = integral_coeff_[i]*integral_[i];
        previous_error[i] = error;
        scalar correction = proportional_correction + differential_correction + integral_correction;
        control_value[i] = min(max(control_value[i] + correction,min_value_[i]),max_value_[i]);
        if (needsUpdate(previous_error[i],i))
        {
            first_iteration_[i] = false;
            offset_value_[i] = 0;
            integral_[i] = integral[i];
            previous_error_[i] = previous_error[i];
            control_value_[i] = control_value[i];
            oldTime_[i] = newTime;
        }
    }
    setNeedsUpdate(previous_error);
    return control_value_;
}

void Foam::feedbackLoopController::setStabilizer(Pair<bool> use_stabilizer,Pair<int> reset_value)
{
    use_stabilizer_ = use_stabilizer;
    positive_feedback_tolerance_ = reset_value;
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
    //max_error_default_ = max_error;
}
void Foam::feedbackLoopController::setMinMaxValue(Pair<scalar> min_value, Pair<scalar> max_value)
{
    forAll(min_value, i)
    {
        if (min_value[i] < max_value[i])
        {
            min_value_[i] = min_value[i];
            max_value_[i] = max_value[i];
        }
    }
}

bool Foam::feedbackLoopController::needsUpdate(scalar error, label i)
{
    return std::abs(error) > std::abs(max_error_[i]*reference_value_[i]);
}

bool Foam::feedbackLoopController::needsUpdate()
{
    return needs_update_;
}

void Foam::feedbackLoopController::setNeedsUpdate(Pair<scalar> error)
{
    needs_update_ = (std::abs(error[0]) > std::abs(max_error_[0]*reference_value_[0]) ||
        std::abs(error[1]) > std::abs(max_error_[1]*reference_value_[1]));
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
