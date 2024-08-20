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
#include <set>
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductingRegionSolvers::conductingRegionSolvers(const Time& runTime)
:
    runTime_(runTime)
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

    //get region characteristic sizes
    characteristicSizes_.setSize(names_.size());
    forAll(names_, i)
    {
        IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                regions_[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        characteristicSizes_[i] =
        (
            physicalProperties.found("Lchar") ?
            dimensionedScalar("Lchar",dimLength,physicalProperties) :
            dimensionedScalar("Lchar",dimLength,0)
        ).value();
    }

    // get write controls / magnetic field update controls
    if (isElectroHarmonic())
    {
        // Maximum allowable magnetic Reynolds number difference comparing
        // to last magnetic field update.
        // This option controls frequency magnetic field updating is called.
        //     (0,inf) - when relative difference in any cell exceeds given value
        //     0     - every iteration
        maxRemDiff_ = readScalar(runTime.controlDict().lookup("maxRemDiff"));

        // Maximum allowable relative field difference in any cell comparing
        // to last magnetic field update.
        // This option controls frequency magnetic field updating is called.
        //     >1  - once
        //     1 - magnetic Reynolds number is used instead
        //     (0,1) - when relative difference in any cell exceeds given value
        //     0     - every iteration
        maxRelDiff_ = readScalar(runTime.controlDict().lookup("maxRelDiff"));
    }
    else
    {
        writeMultiplier_ = readScalar(runTime.controlDict().lookup("writeMultiplier"));
        writeControlDict_ = word(runTime.controlDict().lookup("writeControl"));
        adjustableRunTime_ = (writeControlDict_=="adjustableRunTime");
    }
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
/*
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
*/
//Note: Could make new solver base class with methods for conducting media
//to avoid casting to each of derived class for accessing new common methods.
Foam::solvers::conductingFluid* Foam::conductingRegionSolvers::getFluidPtr_(const word regionName)
{
    if (isFluid(regionName))
    {
        Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
        return dynamic_cast<Foam::solvers::conductingFluid*>(basePtr);
    }
    return nullptr;
}

Foam::solvers::conductingSolid* Foam::conductingRegionSolvers::getSolidPtr_(const word regionName)
{
    if (isSolid(regionName))
    {
        Foam::solver* basePtr = solvers_(regionIdx_[regionName]);
        return dynamic_cast<Foam::solvers::conductingSolid*>(basePtr);
    }
    return nullptr;
}

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
    if (getFluidPtr_(regionName))
    {
        getFluidPtr_(regionName)->solveElectromagnetics();
    }
    else if (getSolidPtr_(regionName))
    {
        getSolidPtr_(regionName)->solveElectromagnetics();
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get solveElectromagnetics!\n"
        << exit(FatalIOError);
    }
}

void Foam::conductingRegionSolvers::electromagneticPredictor(const word regionName)
{
    if (getFluidPtr_(regionName))
    {
        getFluidPtr_(regionName)->electromagneticPredictor();
    }
    else if (getSolidPtr_(regionName))
    {
        getSolidPtr_(regionName)->electromagneticPredictor();
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get electromagneticPredictor!\n"
        << exit(FatalIOError);
    }
}
void Foam::conductingRegionSolvers::calculateGlobalMesh_()
{

    //prepare global mesh

    Info << "Preparing global mesh" << endl;
    int globalCellI = 0;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const cellList cells = regionMesh.cells();
        forAll(cells, cellI)
        {
            //list_c[globalCellI] = cells[cellI];
            localToGlobalID[std::make_pair(regionName,cellI)] = globalCellI++;
        }
    }
    /****************************** */
    /*int points_size = 0;
    int face_size = 0;
    int owner_size = 0;
    int neighbour_size = 0;
    forAll(names_, i)
    {
        fvMesh& regionMesh = regions_[i];
        points_size += regionMesh.points().size();
        face_size += regionMesh.faces().size();
        owner_size += regionMesh.owner().size();
        neighbour_size += regionMesh.neighbour().size();
    }
	const label size_p = points_size;
	const label size_f = face_size;
	const label size_o = owner_size;
	const label size_n = neighbour_size;
    //- list for keeping vector values of field
	List<vector> list_p;
	list_p.setSize(size_p);
	faceList list_f;
	list_f.setSize(size_f);
	labelList list_o;
	list_o.setSize(size_o);
	labelList list_n;
	list_n.setSize(size_n);
    int globalPointI = 0;
    int globalFaceI = 0;
    int globalOwnerI = 0;
    int globalNeighbourI = 0;
    
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
        const labelList owner = regionMesh.owner();
        forAll(owner, cellI)
        {
            list_o[globalOwnerI++] = owner[cellI];
            //localToGlobalID[std::make_pair(regionName,cellI)] = globalCellI++;
        }
        const labelList neighbour = regionMesh.neighbour();
        forAll(neighbour, cellI)
        {
            list_n[globalNeighbourI++] = neighbour[cellI];
            //localToGlobalID[std::make_pair(regionName,cellI)] = globalCellI++;
        }
    }*/

    /****************************** */
	List<vector> list_p;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const faceList faces = regionMesh.faces();
        const pointField points = regionMesh.points();
        forAll(faces, faceI)
        {
            face thisFace = faces[faceI];
            pointField thisFacePoints = thisFace.points(points);
            forAll (thisFacePoints,pointI)
            {
                point thisPoint = thisFacePoints[pointI];
                list_p.append(thisPoint);
            }
        }
    }
    std::map<point,label> pointLabels;
    forAll(list_p, id)
    {
        pointLabels[list_p[id]] = id;
    }
    Info << "Preparing faces " << endl;
	faceList list_f;
    //std::map<label,label> faceLabels;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const faceList faces = regionMesh.faces();
        const pointField points = regionMesh.points();
        forAll(faces, faceI)
        {
            face thisFace = faces[faceI];
            pointField thisFacePoints = thisFace.points(points);
            labelList thisFaceLabels;
            forAll (thisFacePoints,pointI)
            {
                point thisPoint = thisFacePoints[pointI];
                thisFaceLabels.append(pointLabels[thisPoint]);
            }
            face newFace(thisFaceLabels);
            //faceLabels[faceI] = list_f.size();
            list_f.append(newFace);
        }
    }

    /*Pout << "Preparing cells " << endl;
    cellList list_c;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const faceList faces = regionMesh.faces();
        const pointField points = regionMesh.points();
        const cellList cells = regionMesh.cells();
        forAll(cells, cellI)
        {
            cell thisCell = cells[cellI];
            labelList thisCellLabels;
            forAll (thisCell,faceI)
            {
                thisCellLabels.append(faceLabels[thisCell[faceI]]);
            }
            cell newCell(thisCellLabels);
            list_c.append(newCell);
        }
    }*/
    labelList list_o;
    labelList list_n;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const labelList owner = regionMesh.owner();
        //const cellList cells = regionMesh.cells();
        forAll(owner, cellI)
        {
            list_o.append(localToGlobalID[std::make_pair(regionName,owner[cellI])]);
        }
        const labelList neighbour = regionMesh.neighbour();
        forAll(neighbour, cellI)
        {
            list_n.append(localToGlobalID[std::make_pair(regionName,neighbour[cellI])]);
        }
    }
    /**************************
    Needs ordering
    ****************/

    //Pout << "points " << endl;
    //Pout << list_p << endl;
    //Pout << "faces " << endl;
    //Pout << list_f << endl;
    //Pout << "cells " << endl;
    //Pout << list_c << endl;
    //- reference list used for setting values to field  
	//const vectorUList& Ulist_v(list_v);
    //field.internalField() = vectorField(Ulist_v);


    Pout << "assigning mesh " << endl;
    globalMesh_ = new 
    fvMesh
        (
            IOobject
            (
                "dummyRegion",//fvMesh::defaultRegion,//
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ
            ),
            pointField(list_p),
            faceList(list_f),
            labelList(list_o),
            labelList(list_n),
            //cellList(list_c),
            false
        );
    Pout << "global mesh " << endl;
    const Foam::fvMesh& globalMesh = globalMesh_;
    //Pout << "global write " << endl;
    //globalMesh.write();
    /*
    Pout << "localToGlobalID " << endl;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        forAll(regionMesh.cells(), cellI)
        {
            Pout << cellI << "   "
            << localToGlobalID[std::make_pair(regionName,cellI)] << endl;
        }
    }
    const cellList cells = glob.cells();
    forAll(cells, cellI)
    {
        Pout << cellI << "  " << cells[cellI] << endl;
    }*/
    point tmpPoint = globalMesh.points()[0];
    Pout << "point " << endl;
    Pout << "point test: " <<   globalMesh.findCell(tmpPoint) << endl;


    if (!globalMesh.objectRegistry::foundObject<IOdictionary>("fvSchemes"))
    {
        IOdictionary* fPtr
        (
            new IOdictionary
            (
                
            IOobject
            (
                "fvSchemes",
                globalMesh.time().system(),
                globalMesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

}
const Foam::fvMesh& Foam::conductingRegionSolvers::globalMesh()
{
    if (globalMesh_.empty())
    {
        calculateGlobalMesh_();
    }

    if (globalMesh_.valid())
    {
        return globalMesh_;
    }
    else
    {
        FatalIOError
        << "Failed to get global mesh!\n"
        << exit(FatalIOError);
    }
}
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolvers::setJ(volVectorField& globalField, bool imaginary)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getFluidPtr_(regionName))
        {
            volVectorField& regionField = getFluidPtr_(regionName)->getJ(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
        }
        else if (getSolidPtr_(regionName))
        {
            volVectorField& regionField = getSolidPtr_(regionName)->getJ(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
        }
        else
        {
            FatalIOError
            << " region " << regionName << " solver is not " << fluidSolverName_
            << " or " << solidSolverName_ << "!\n" << "Cannot getJ!\n"
            << exit(FatalIOError);
        }
    }
}
//Assigns fluid and solid region values from global to each region field
void Foam::conductingRegionSolvers::setB(volVectorField& globalField, bool imaginary)
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getFluidPtr_(regionName))
        {
            volVectorField& regionField = getFluidPtr_(regionName)->getB(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
        }
        else if (getSolidPtr_(regionName))
        {
            volVectorField& regionField = getSolidPtr_(regionName)->getB(imaginary);
            vectorGlobalToField_(globalField,regionField,regionName);
        }
        else
        {
            FatalIOError
            << " region " << regionName << " solver is not " << fluidSolverName_
            << " or " << solidSolverName_ << "!\n" << "Cannot getB!\n"
            << exit(FatalIOError);
        }
    }
}
//Assigns U_old = U
void Foam::conductingRegionSolvers::storeU(const word regionName)
{
    if (getFluidPtr_(regionName))
    {
        getFluidPtr_(regionName)->storeU();
    }
}
//returns read-only access to electro module
const Foam::electromagneticModel& Foam::conductingRegionSolvers::getElectro(const word regionName)
{
    if (isFluid(regionName))
    {
        return getFluid(regionName).electro;
    }
    else if (isSolid(regionName))
    {
        return getSolid(regionName).electro;
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get electromagneticPredictor!\n"
        << exit(FatalIOError);
    }
}
const Foam::solvers::conductingFluid& Foam::conductingRegionSolvers::getFluid(const word regionName)
{
    if (getFluidPtr_(regionName))
    {
        const Foam::solvers::conductingFluid& fluidRef(*getFluidPtr_(regionName));
        return fluidRef;
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get fluid!\n"
        << exit(FatalIOError);
    }
}
const Foam::solvers::conductingSolid& Foam::conductingRegionSolvers::getSolid(const word regionName)
{
    if (getSolidPtr_(regionName))
    {
        const Foam::solvers::conductingSolid& solidRef(*getSolidPtr_(regionName));
        return solidRef;
    }
    else
    {
        FatalIOError
        << " region " << regionName << " solver is not " << fluidSolverName_
        << " or " << solidSolverName_ << "!\n" << "Cannot get solid!\n"
        << exit(FatalIOError);
    }
}

//Assigns single fluid/solid region values from region to global field
void Foam::conductingRegionSolvers::scalarFieldToGlobal(volScalarField& global,const volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[localToGlobalID[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}

void Foam::conductingRegionSolvers::vectorFieldToGlobal(volVectorField& global,const volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[localToGlobalID[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}
bool Foam::conductingRegionSolvers::isFluid(const word regionName)
{
    return names_[regionIdx_[regionName]].second() == fluidSolverName_;
}

bool Foam::conductingRegionSolvers::isSolid(const word regionName)
{
    return names_[regionIdx_[regionName]].second() == solidSolverName_;
}

void Foam::conductingRegionSolvers::setCorrectElectromagnetics()
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (getFluidPtr_(regionName))
        {
            getFluidPtr_(regionName)->setCorrectElectromagnetics();
        }
        else if (getSolidPtr_(regionName))
        {
            getSolidPtr_(regionName)->setCorrectElectromagnetics();
        }
        else
        {
            FatalIOError
            << " region " << regionName << " solver is not " << fluidSolverName_
            << " or " << solidSolverName_ << "!\n" << "Cannot get electromagneticPredictor!\n"
            << exit(FatalIOError);
        }
    }
}

bool Foam::conductingRegionSolvers::isElectroHarmonic()
{
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        if (isFluid(regionName) || isSolid(regionName))
        {
            return getElectro(regionName).isComplex();
        }
    }
    return false;
}

bool Foam::conductingRegionSolvers::updateMagneticField()
{
    bool doUpdate = false;
    if (isElectroHarmonic())
    {
        scalar maxRemDiff_local = SMALL;        
        scalar maxRelDiff_local = SMALL;

        forAll(names_, i)
        {
            const word& regionName = names_[i].first();
            if (isFluid(regionName))
            {
                const volVectorField& U = getFluid(regionName).U;
                const volVectorField& U_old = getFluid(regionName).U_old;
                maxRemDiff_local = max(
                    mu_0 * characteristicSizes_[i] *
                    max(getElectro(regionName).sigmaInv()*mag(U_old-U)).value(),
                    maxRemDiff_local);        

                maxRelDiff_local = max(
                    (max(mag(U_old-U)/(average(mag(U))+smallU))).value(),
                    maxRemDiff_local);
            }
        }

        if((maxRelDiff_local>maxRelDiff_ || maxRelDiff_<SMALL) && maxRelDiff_+SMALL<=1.0) {
            doUpdate = true;
        }
        else if(maxRemDiff_local>maxRemDiff_ && maxRelDiff_-SMALL<=1.0) {
            doUpdate = true;
        }
    }
    else
    {
        if (adjustableRunTime_)
        {
            if (runTime_.writeTime()) doUpdate = true;
        }
        else
        {
            writeCounter_++;
            if ( (writeCounter_ % writeMultiplier_) == 0 && runTime_.run()) doUpdate = true;
        }
    }
    return doUpdate;
}

void Foam::conductingRegionSolvers::calcTemperatureGradient(const word regionName)
{
    Pout << "calc gradient" << endl;
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime_.constant(),
            mesh(regionName),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    if
    ( 
        physicalProperties.found("temperature_multiplier") &&
        physicalProperties.found("temperature_addition")
    )
    {
        dimensionedVector temperature_multiplier
        (
            "temperature_multiplier",
            dimTemperature/dimLength,
            physicalProperties
        );
        dimensionedScalar temperature_addition
        (
            "temperature_addition",
            dimTemperature,
            physicalProperties
        );

        if (getSolidPtr_(regionName))
        {
            volScalarField& T = getSolidPtr_(regionName)->getTemperature();
            T = (mesh(regionName).C() & temperature_multiplier) +  temperature_addition;
            T.write();
        }
        else if (getFluidPtr_(regionName))
        {
            volScalarField& T = getFluidPtr_(regionName)->getTemperature();
            T = (mesh(regionName).C() & temperature_multiplier) +  temperature_addition;
            T.write();
        }
        /**/
        else
        {
            Info << "Warning: region " << regionName << " solver is not " << fluidSolverName_
            << " or " << solidSolverName_ << "!\n" << "Cannot set temperature gradient!\n"; 
        }
        /**/
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

Foam::solver& Foam::conductingRegionSolvers::operator()(const word regionName)
{
    setPrefix(regionIdx_[regionName]);
    return solvers_[regionIdx_[regionName]];
}
/*
Foam::solver* Foam::conductingRegionSolvers::operator()(const label i)
{
    return solvers_(i);
}*/


// ************************************************************************* //
