#include "createFluidFields.H"
#include "createSolidFields.H"
Info<< "Initializing field UGlobal\n" << endl;
volVectorField UGlobal
(
    IOobject
    (
        "U",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::vector(0,0,0))
);
//Use only fluid mesh if single fluid region
//Use global mesh if multiple fluid regions
volVectorField Usent = fluidRegions.size() == 1 ? UFluid[0] : UGlobal;

volVectorField JGlobal
(
    IOobject
    (
        "J",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, -1, -3, 0, 0, 0, 0), Foam::vector(0,0,0))
);
Info<< "Initializing field BGlobal\n" << endl;
volVectorField BGlobal
(
    IOobject
    (
        "B",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, -1, -3, 0, 0, 0, 0), Foam::vector(0,0,0))
);
