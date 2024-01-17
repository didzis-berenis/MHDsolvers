regionProperties rp(runTime);

#include "createFluidMeshes.H"
#include "createSolidMeshes.H"

Foam::fvMesh meshGlobal
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
