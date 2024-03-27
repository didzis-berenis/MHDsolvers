
Info<< "Initializing field J\n" << endl;
volVectorField J
(
    IOobject
    (
        "J",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Electric potential field
Info<< "Reading field PotE\n" << endl;
volScalarField PotE
(
    IOobject
    (
        "PotE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Electric potential reference point
label PotERefCell = 0;
scalar PotERefValue = 0.0;
setRefCell(PotE, solver.pimple.dict(), PotERefCell, PotERefValue);
mesh.schemes().setFluxRequired(PotE.name());

int writeMultiplier( readScalar(runTime.controlDict().lookup("writeMultiplier")) );
long long int writeCounter = 0;
word writeControlDict(runTime.controlDict().lookup("writeControl"));
const bool adjustableRunTime = (writeControlDict=="adjustableRunTime");
