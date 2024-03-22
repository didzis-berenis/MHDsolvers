//regionProperties rp(runTime);

//#include "createFluidMeshes.H"
//#include "createSolidMeshes.H"
// Read global mesh of all regions
fvMesh meshGlobal
(
    IOobject
    (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);

IOList<label> cellToRegion
(
    IOobject
    (
        "cellToRegion",
        runTime.constant(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )

);

int offset = 0;
forAll(meshGlobal.cellZones(), zoneI)
{
    const cellZone& cZone = meshGlobal.cellZones()[zoneI];
    regionOffset[cZone.name()] = offset;
    offset += cZone.size();
    regionID[cZone.name()] = zoneI;
}

forAll(regionSolverNames, i)
{
    const word& solverName = regionSolverNames[i].second();
    if (solverName == "conductingFluid")
    {
        fluidID[solverName] = i;
    }
    if (solverName == "conductingSolid")
    {
        solidID[solverName] = i;
    }
}

