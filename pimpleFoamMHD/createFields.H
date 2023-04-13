#include "createRDeltaT.H"

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

#include "createPhi.H"

pressureReference pressureReference(p, pimple.dict());

mesh.schemes().setFluxRequired(p.name());


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

// Maximum allowable relative field difference in any cell comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     1     - once
//     (0,1) - when relative difference in any cell exceeds given value
//     0     - every iterationInfo<< "Reading physicalProperties\n" << endl;
scalar maxRelDiff(readScalar(physicalProperties.lookup("maxRelDiff")));

