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

#include "conductingRegionSolvers.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::conductingRegionSolvers(const Time& runTime)
{

    if (runTime.controlDict().found("regionSolvers"))
    {
        const dictionary& conductingRegionSolversDict =
            runTime.controlDict().subDict("regionSolvers");

        forAllConstIter(dictionary, conductingRegionSolversDict, iter)
        {
            const word regionName(iter().keyword());
            const word solverName(iter().stream());

            names_.append(Pair<word>(regionName, solverName));
        }
    }
    else
    {
        FatalIOErrorInFunction(runTime.controlDict())
            << "regionSolvers list missing from "
            << runTime.controlDict().name()
            << exit(FatalIOError);
    }

    regions_.setSize(names_.size());
    solvers_.setSize(names_.size());
    prefixes_.setSize(names_.size());

    string::size_type nRegionNameChars = 0;

    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        const word& solverName = names_[i].second();
        regionIdx_[regionName] = i;

        // Load the solver library
        solver::load(solverName);

        regions_.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        //Foam::autoPtr<Foam::solver> solverPtr = solver::New(solverName, regions_[i]);
        solvers_.set(i, solver::New(solverName, regions_[i]));

        prefixes_[i] = regionName;
        nRegionNameChars = max(nRegionNameChars, regionName.size());
    }

    nRegionNameChars++;

    prefix0_.append(nRegionNameChars, ' ');

    forAll(names_, i)
    {
        prefixes_[i].append(nRegionNameChars - prefixes_[i].size(), ' ');
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::~conductingRegionSolvers()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conductingRegionSolvers::setGlobalPrefix() const
{
    Sout.prefix() = prefix0_;
}


void Foam::conductingRegionSolvers::setPrefix(const label i) const
{
    Sout.prefix() = prefixes_[i];
}


void Foam::conductingRegionSolvers::resetPrefix() const
{
    Sout.prefix() = string::null;
}


Foam::fvMesh& Foam::conductingRegionSolvers::mesh(const word regionName)
{
    return regions_[regionIdx_[regionName]];
}

Foam::solvers::conductingFluid* Foam::conductingRegionSolvers::getFluidPtr_(const word regionName)
{
    const word& solverName = names_[regionIdx_[regionName]].second();
    Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
    if (solverName == fluidSolverName_)
    {
        return dynamic_cast<Foam::solvers::conductingFluid*>(basePtr);
    }
    return nullptr;
}

Foam::solvers::conductingSolid* Foam::conductingRegionSolvers::getSolidPtr_(const word regionName)
{
    const word& solverName = names_[regionIdx_[regionName]].second();
    Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
    if (solverName == solidSolverName_)
    {
        return dynamic_cast<Foam::solvers::conductingSolid*>(basePtr);
    }
    return nullptr;
}

// Return region fvMesh
void Foam::conductingRegionSolvers::setJxB(const word regionName, Foam::volVectorField& field)
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_(regionName);
    if (fluidPtr)
    {
        fluidPtr->setJxB(field);
    }
    else
    {
        Info << "Warning: region " << regionName << " solver is not " << fluidSolverName_
        << "!\n" << "Cannot set JxB field!\n"; 
    }
}

// Return region fvMesh
void Foam::conductingRegionSolvers::setJJsigma(const word regionName, Foam::volScalarField& field)
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_(regionName);
    Foam::solvers::conductingSolid* solidPtr = getSolidPtr_(regionName);
    if (fluidPtr)
    {
        fluidPtr->setJJsigma(field);
    }
    else if (solidPtr)
    {
        solidPtr->setJJsigma(field);
    }
    else
    {
        Info << "Warning: region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot set JJsigma field!\n"; 
    }
}
/*
// Return region fvMesh
Foam::volVectorField& Foam::conductingRegionSolvers::getVelocity(const word regionName)
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_(regionName);
    if (fluidPtr)
    {
        return fluidPtr->getVelocity();
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << "!\n" << "Cannot get Velocity field!\n"
        << exit(FatalIOError);
    }
}

// Return region fvMesh
Foam::volScalarField& Foam::conductingRegionSolvers::getPressure(const word regionName)
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_(regionName);
    if (fluidPtr)
    {
        return fluidPtr->getPressure();
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << "!\n" << "Cannot get Pressure field!\n"
        << exit(FatalIOError);
    }
}

// Return region fvMesh
Foam::volScalarField& Foam::conductingRegionSolvers::getTemperature(const word regionName)
{
    Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_(regionName);
    Foam::solvers::conductingSolid* solidPtr = getSolidPtr_(regionName);
    if (fluidPtr)
    {
        return fluidPtr->getTemperature();
    }
    else if (solidPtr)
    {
        return solidPtr->getTemperature();
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Temperature field!\n"
        << exit(FatalIOError);
    }
}
*/
bool Foam::conductingRegionSolvers::isFluid(const word regionName)
{
    if (names_[regionIdx_[regionName]].second() == fluidSolverName_)
    {
        return true;
    }
    return false;
}

bool Foam::conductingRegionSolvers::isSolid(const word regionName)
{
    if (names_[regionIdx_[regionName]].second() == solidSolverName_)
    {
        return true;
    }
    return false;
}

void Foam::conductingRegionSolvers::setCorrectElectromagnetics()
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        Foam::solvers::conductingFluid* fluidPtr = getFluidPtr_(regionName);
        Foam::solvers::conductingSolid* solidPtr = getSolidPtr_(regionName);
        if (fluidPtr)
        {
            fluidPtr->setCorrectElectromagnetics();
        }
        else if (solidPtr)
        {
            solidPtr->setCorrectElectromagnetics();
        }
    }
}

Foam::List<Foam::Pair<Foam::word>> Foam::conductingRegionSolvers::getNames()
{
    return names_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::solver& Foam::conductingRegionSolvers::operator[](const label i)
{
    setPrefix(i);
    return solvers_[i];
}
/*
Foam::solver* Foam::conductingRegionSolvers::operator()(const label i)
{
    return solvers_(i);
}*/


// ************************************************************************* //
