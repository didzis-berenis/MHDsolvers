
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
