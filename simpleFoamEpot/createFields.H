Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field U_old\n" << endl;
volVectorField U_old
(
    IOobject
    (
        "U_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);
Info<< "Initializing field JxB\n" << endl;
volVectorField JxB
(
    IOobject
    (
        "JxB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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


#include "createPhi.H"

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

pressureReference pressureReference(p, simple.dict());
mesh.schemes().setFluxRequired(p.name());

//Electric potential reference point
label PotEreRefCell = 0;
scalar PotEreRefValue = 0.0;
setRefCell(PotEre, simple.dict(), PotEreRefCell, PotEreRefValue);
mesh.schemes().setFluxRequired(PotEre.name());

label PotEimRefCell = 0;
scalar PotEimRefValue = 0.0;
setRefCell(PotEim, simple.dict(), PotEimRefCell, PotEimRefValue);
mesh.schemes().setFluxRequired(PotEim.name());


autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
dimensionedScalar rho("rho", dimensionSet(1, -3, 0, 0, 0, 0, 0), physicalProperties);

//Electrical conductivity
dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    physicalProperties
);

// Maximum allowable magnetic Reynolds number difference comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     (0,inf) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRemDiff(readScalar(physicalProperties.lookup("maxRemDiff")));

//
// Maximum allowable relative field difference in any cell comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     >1  - once
//     1 - magnetic Reynolds number is used instead
//     (0,1) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRelDiff(readScalar(physicalProperties.lookup("maxRelDiff")));

