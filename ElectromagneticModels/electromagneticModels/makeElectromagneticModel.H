/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

InClass
    Foam::solidThermo

Description
    Macros for creating solid thermo packages

\*---------------------------------------------------------------------------*/

#include "transientElectromagneticModel.H"
#include "harmonicElectromagneticModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeElectromagneticModels()                                                     \
                                                                                        \
    namespace Foam                                                                      \
    {                                                                                   \
                                                                                        \
        defineNamedTemplateTypeNameAndDebug(harmonicElectromagneticModel, 0);           \
                                                                                        \
        defineTemplateRunTimeSelectionTable(harmonicElectromagneticModel, dictionary);  \
                                                                                        \
        addToRunTimeSelectionTable                                                      \
        (                                                                               \
            harmonicElectromagneticModel,                                               \
            dictionary                                                                  \
        );                                                                              \
                                                                                        \
                                                                                        \
        defineNamedTemplateTypeNameAndDebug(transientElectromagneticModel, 0);          \
                                                                                        \
        defineTemplateRunTimeSelectionTable(transientElectromagneticModel, dictionary); \
                                                                                        \
        addToRunTimeSelectionTable                                                      \
        (                                                                               \
            transientElectromagneticModel,                                              \
            dictionary                                                                  \
        );                                                                              \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
