
Info<< "Initializing field JreGlobal\n" << endl;
volVectorField JreGlobal
(
    IOobject
    (
        "Jre",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, -1, -3, 0, 0, 0, 0), Foam::vector(0,0,0))
);
Info<< "Initializing field BreGlobal\n" << endl;
volVectorField BreGlobal
(
    IOobject
    (
        "Bre",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, -1, -3, 0, 0, 0, 0), Foam::vector(0,0,0))
);
Info<< "Initializing field JimGlobal\n" << endl;
volVectorField JimGlobal
(
    IOobject
    (
        "Jim",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, -1, -3, 0, 0, 0, 0), Foam::vector(0,0,0))
);
Info<< "Initializing field BimGlobal\n" << endl;
volVectorField BimGlobal
(
    IOobject
    (
        "Bim",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, -1, -3, 0, 0, 0, 0), Foam::vector(0,0,0))
);
