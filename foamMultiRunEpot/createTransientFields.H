
volVectorField JGlobal
(
    IOobject
    (
        "J",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field BGlobal\n" << endl;
volVectorField BGlobal
(
    IOobject
    (
        "B",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field PotEGlobal\n" << endl;
volScalarField PotEGlobal
(
    IOobject
    (
        "PotE",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
meshGlobal.schemes().setFluxRequired(PotEGlobal.name());
