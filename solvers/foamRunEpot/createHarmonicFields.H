
Info<< "Initializing field Jre\n" << endl;
volVectorField Jre
(
    IOobject
    (
        "Jre",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field Jim\n" << endl;
volVectorField Jim
(
    IOobject
    (
        "Jim",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field Bre\n" << endl;
volVectorField Bre
(
    IOobject
    (
        "Bre",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field Bim\n" << endl;
volVectorField Bim
(
    IOobject
    (
        "Bim",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Electric potential field
Info<< "Reading field PotEre\n" << endl;
volScalarField PotEre
(
    IOobject
    (
        "PotEre",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field PotEim\n" << endl;
volScalarField PotEim
(
    IOobject
    (
        "PotEim",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Electric potential reference point
label PotEreRefCell = 0;
scalar PotEreRefValue = 0.0;
setRefCell(PotEre, solver.pimple.dict(), PotEreRefCell, PotEreRefValue);
mesh.schemes().setFluxRequired(PotEre.name());

label PotEimRefCell = 0;
scalar PotEimRefValue = 0.0;
setRefCell(PotEim, solver.pimple.dict(), PotEimRefCell, PotEimRefValue);
mesh.schemes().setFluxRequired(PotEim.name());

//Characteristic length
dimensionedScalar Lchar
(
    "Lchar",
    dimensionSet(0, 1, 0, 0, 0, 0, 0),
    physicalProperties
);

double Rem0 = 4*3.14159*(std::pow(10,-7))*sigma.value()*Lchar.value();
Info<< "Rem0 = " << Rem0 << endl;

// Maximum allowable magnetic Reynolds number difference comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     (0,inf) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRemDiff(readScalar(runTime.controlDict().lookup("maxRemDiff")));

//
// Maximum allowable relative field difference in any cell comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     >1  - once
//     1 - magnetic Reynolds number is used instead
//     (0,1) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRelDiff(readScalar(runTime.controlDict().lookup("maxRelDiff")));
