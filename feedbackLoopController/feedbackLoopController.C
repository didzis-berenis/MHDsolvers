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
    //fvMesh& mesh,
    const scalar& proportional_coeff,
    const scalar& differential_coeff,
    const scalar& integral_coeff,
    const scalar& required_output
)
:
    proportional_coeff_(proportional_coeff),
    differential_coeff_(differential_coeff),
    integral_coeff_(integral_coeff),
    required_output_(required_output),
    previous_error_(0),
    integral_(0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::feedbackLoopController::~feedbackLoopController()
{}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::scalar Foam::feedbackLoopController::calculateCorrection(scalar present_output, scalar deltaT)
{
    scalar error = required_output_-present_output;
    scalar proportional_correction = proportional_coeff_*error;
    scalar differential_correction = 0;
    scalar integral_correction = 0;
    if (deltaT > 0)
    {
        differential_correction = differential_coeff_*(error-previous_error_)/deltaT;
        integral_ += error*deltaT;
        integral_correction = integral_coeff_*integral_;
    }
    return proportional_correction + differential_correction + integral_correction;
}

// ************************************************************************* //
