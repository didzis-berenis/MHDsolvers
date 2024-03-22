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

#include "regionSolvers.H"
#include "solver.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::regionSolvers(const Time& runTime)
{

    if (runTime.controlDict().found("regionSolvers"))
    {
        const dictionary& regionSolversDict =
            runTime.controlDict().subDict("regionSolvers");

        forAllConstIter(dictionary, regionSolversDict, iter)
        {
            const word regionName(iter().keyword());
            const word solverName(iter().stream());

            names_.append(Pair<word>(regionName, solverName));
        }
    }
    else
    {
        // Partial backward-compatibility
        // Converts the regions entry in the regionProperties dictionary into
        // the regionSolvers list
        // Only supports fluid and solid regions

        typeIOobject<IOdictionary> regionPropertiesHeader
        (
            IOobject
            (
                "regionProperties",
                runTime.time().constant(),
                runTime.db(),
                IOobject::MUST_READ
            )
        );

        if (regionPropertiesHeader.headerOk())
        {
            HashTable<wordList> regions
            (
                IOdictionary(regionPropertiesHeader).lookup("regions")
            );

            if (regions.found("solid"))
            {
                const wordList& fluidRegions = regions["solid"];
                forAll(fluidRegions, i)
                {
                    names_.append
                    (
                        Pair<word>(fluidRegions[i], "conductingSolid")
                    );
                }
            }

            if (regions.found("fluid"))
            {
                const wordList& fluidRegions = regions["fluid"];
                forAll(fluidRegions, i)
                {
                    names_.append
                    (
                        Pair<word>(fluidRegions[i], "conductingFluid")
                    );
                }
            }
        }
        else
        {
            FatalIOErrorInFunction(runTime.controlDict())
                << "regionSolvers list missing from "
                << runTime.controlDict().name()
                << exit(FatalIOError);
        }
    }

    regions_.setSize(names_.size());
    solvers_.setSize(names_.size());
    prefixes_.setSize(names_.size());

    string::size_type nRegionNameChars = 0;

    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        const word& solverName = names_[i].second();
        //int fluidCounter = 0;
        //int solidCounter = 0;

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
        solvers_.set(i, solver::New(solverName, regions_[i]));/*
        if (solverName == "conductingFluid")
        {
            Pout<< nl << "Setting fluid\n" << endl;
            //Foam::autoPtr<Foam::solvers::conductingFluid> electromagneticSolverPtr(dynamic_cast<Foam::solvers::conductingFluid*>(solverPtr.ptr()));
            //fluids_.set(fluidCounter,electromagneticSolverPtr);
            fluidCounter++;
        }
        if (solverName == "conductingSolid")
        {
            Pout<< nl << "Setting solid\n" << endl;
            //Foam::autoPtr<Foam::solvers::conductingSolid> electromagneticSolverPtr(dynamic_cast<Foam::solvers::conductingSolid*>(solverPtr.ptr()));
            //solids_.set(solidCounter,electromagneticSolverPtr);
            solidCounter++;
        }*/

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

Foam::regionSolvers::~regionSolvers()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::setGlobalPrefix() const
{
    Sout.prefix() = prefix0_;
}


void Foam::regionSolvers::setPrefix(const label i) const
{
    Sout.prefix() = prefixes_[i];
}


void Foam::regionSolvers::resetPrefix() const
{
    Sout.prefix() = string::null;
}


Foam::fvMesh& Foam::regionSolvers::mesh(const label i)
{
    return regions_[i];
}

// Return region fvMesh
Foam::solvers::conductingFluid* Foam::regionSolvers::fluids(const label i)
{
    return dynamic_cast<Foam::solvers::conductingFluid*>(solvers_(i));
}

// Return region fvMesh
Foam::solvers::conductingSolid* Foam::regionSolvers::solids(const label i)
{
    return dynamic_cast<Foam::solvers::conductingSolid*>(solvers_(i));
}


Foam::List<Foam::Pair<Foam::word>> Foam::regionSolvers::getNames()
{
    return names_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::solver& Foam::regionSolvers::operator[](const label i)
{
    setPrefix(i);
    return solvers_[i];
}
/*
Foam::solver* Foam::regionSolvers::operator()(const label i)
{
    return solvers_(i);
}*/


// ************************************************************************* //
