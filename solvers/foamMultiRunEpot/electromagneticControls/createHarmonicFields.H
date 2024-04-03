
Info<< "Initializing field JreGlobal\n" << endl;
volVectorField JreGlobal
(
    IOobject
    (
        "Jre",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field JimGlobal\n" << endl;
volVectorField JimGlobal
(
    IOobject
    (
        "Jim",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field BreGlobal\n" << endl;
volVectorField BreGlobal
(
    IOobject
    (
        "Bre",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field BimGlobal\n" << endl;
volVectorField BimGlobal
(
    IOobject
    (
        "Bim",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field PotEreGlobal\n" << endl;
volScalarField PotEreGlobal
(
    IOobject
    (
        "PotEre",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field PotEimGlobal\n" << endl;
volScalarField PotEimGlobal
(
    IOobject
    (
        "PotEim",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal
);

meshGlobal.schemes().setFluxRequired(PotEreGlobal.name());

meshGlobal.schemes().setFluxRequired(PotEimGlobal.name());
