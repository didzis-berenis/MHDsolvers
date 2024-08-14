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

    int points_size = 0;
    int face_size = 0;
    int cell_size = 0;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        points_size += regionMesh.points().size();
        face_size += regionMesh.faces().size();
        cell_size += regionMesh.cells().size();
    }
	const label size_p = points_size;	
	const label size_f = face_size;	
	const label size_c = cell_size;	
    //- list for keeping vector values of field
	List<vector> list_p;
	list_p.setSize(size_p);
	faceList list_f;
	list_f.setSize(size_f);
	cellList list_c;
	list_c.setSize(size_c);
    int globalPointI = 0;
    int globalFaceI = 0;
    int globalCellI = 0;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const pointField points = regionMesh.points();
        forAll(points, pointI)
        {
            list_p[globalPointI++] = points[pointI];
        }
        const faceList faces = regionMesh.faces();
        forAll(faces, faceI)
        {
            list_f[globalFaceI++] = faces[faceI];
        }
        const cellList cells = regionMesh.cells();
        forAll(cells, cellI)
        {
            list_c[globalCellI] = cells[cellI];
            localToGlobalID[std::make_pair(regionName,cellI)] = globalCellI++;
        }
    }
    //- reference list used for setting values to field  
	//const vectorUList& Ulist_v(list_v);
    //field.internalField() = vectorField(Ulist_v);

    globalMesh_ = new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            pointField(list_p),
            faceList(list_f),
            cellList(list_c),
            false
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::~conductingRegionSolvers()
{}

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

//Assigns single fluid/solid region values from global field to region
void Foam::conductingRegionSolvers::scalarGlobalToField_(volScalarField& global,volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[localToGlobalID[std::make_pair(regionName,cellI)]];
    }
}

void Foam::conductingRegionSolvers::vectorGlobalToField_(volVectorField& global,volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[localToGlobalID[std::make_pair(regionName,cellI)]];
    }
}

Foam::electroBase* Foam::conductingRegionSolvers::getElectroBasePtr_(const word regionName)
{
    const word& solverName = names_[regionIdx_[regionName]].second();
    Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
    if (solverName == fluidSolverName_ || solverName == solidSolverName_)
    {
        return dynamic_cast<Foam::electroBase*>(basePtr);
    }
    return nullptr;
}
/*
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
*/

// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

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

void Foam::conductingRegionSolvers::solveElectromagnetics(const word regionName)
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_(regionName);
    if (electroBasePtr)
    {
        electroBasePtr->solveElectromagnetics();
    }
    /**/
    else
    {
        Info << "Warning: region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot solve electromagnetics!\n"; 
    }
    /**/
}

const Foam::fvMesh& Foam::conductingRegionSolvers::globalMesh()
{
    return globalMesh_;
}

/*
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

//Assigns fluid and solid region values from each region to global field
/*
void Foam::conductingRegionSolvers::vectorMultiRegionToGlobal(volVectorField& global)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            global[localToGlobalID[std::make_pair(regionNames[i],cellI)]] = regions[i][cellI];
        }
    }
}

void Foam::conductingRegionSolvers::scalarMultiRegionToGlobal(volScalarField& global)
{
    forAll(names_, i)
    {
        forAll(regions[i], cellI)
        {
            global[localToGlobalID[std::make_pair(regionNames[i],cellI)]] = regions[i][cellI];
        }
    }
}
*/
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolvers::setJ(volVectorField& globalField, bool imaginary)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        Foam::electroBase* electroBasePtr = getElectroBasePtr_(regionName);
        if (electroBasePtr)
        {
            volVectorField& regionField = electroBasePtr->getJ(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
        }
    }
}
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolvers::setB(volVectorField& globalField, bool imaginary)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        Foam::electroBase* electroBasePtr = getElectroBasePtr_(regionName);
        if (electroBasePtr)
        {
            volVectorField& regionField = electroBasePtr->getB(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
        }
    }
}
//returns read-only access to electro module
const Foam::electromagneticModel& Foam::conductingRegionSolvers::getElectro(const word regionName)
{
    Foam::electroBase* electroBasePtr = getElectroBasePtr_(regionName);
    if (electroBasePtr)
    {
        return electroBasePtr->electro;
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get Electric module!\n"
        << exit(FatalIOError);
    }
}
//Assigns single fluid/solid region values from region to global field
void Foam::conductingRegionSolvers::scalarFieldToGlobal(volScalarField& global,volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[localToGlobalID[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}

void Foam::conductingRegionSolvers::vectorFieldToGlobal(volVectorField& global,volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[localToGlobalID[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}
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
        Foam::electroBase* electroBasePtr = getElectroBasePtr_(regionName);
        if (electroBasePtr)
        {
            electroBasePtr->setCorrectElectromagnetics();
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
