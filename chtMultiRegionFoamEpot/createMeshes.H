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
