regionProperties rp(runTime);

#include "createFluidMeshes.H"
#include "createSolidMeshes.H"
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

std::map<word, label> regionOffset;
std::map<word, label> regionID;

int offset = 0;
forAll(meshGlobal.cellZones(), zoneI)
{
    const cellZone& cZone = meshGlobal.cellZones()[zoneI];
    regionOffset[cZone.name()] = offset;
    offset += cZone.size();
    regionID[cZone.name()] = zoneI;
}

// globalToRegionIDFluid and globalToRegionIDSolid could be used
// to map fluid/solid region id from global IDs, which are stored in cellToRegion

//std::map<label, label> globalToRegionIDFluid;
//std::map<label, label> globalToRegionIDSolid;

forAll(fluidNames, i)
{ 
    fluidOffset[i] = regionOffset[fluidNames[i]];
    //globalToRegionIDFluid.insert(std::pair<label, label>(regionID[fluidNames[i]], i));
}
forAll(solidNames, i)
{ 
    solidOffset[i] = regionOffset[solidNames[i]];
    //globalToRegionIDSolid.insert(std::pair<label, label>(regionID[fluidNames[i]], i));
}

