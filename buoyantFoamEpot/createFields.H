#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
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


Info<< "Initializing field JJsigma\n" << endl;
volScalarField JJsigma
(
    IOobject
    (
        "JJsigma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#if (ELMER_TIME == HARMONIC)
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

#elif (ELMER_TIME == TRANSIENT)
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
#endif

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    thermo.incompressible()
);

mesh.schemes().setFluxRequired(p_rgh.name());

//Electric potential reference point
label PotEreRefCell = 0;
scalar PotEreRefValue = 0.0;
setRefCell(PotEre, pimple.dict(), PotEreRefCell, PotEreRefValue);
mesh.schemes().setFluxRequired(PotEre.name());

label PotEimRefCell = 0;
scalar PotEimRefValue = 0.0;
setRefCell(PotEim, pimple.dict(), PotEimRefCell, PotEimRefValue);
mesh.schemes().setFluxRequired(PotEim.name());

hydrostaticInitialisation
(
    p_rgh,
    p,
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
#include "checkRadiationModel.H"

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

//Electrical conductivity
dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    physicalProperties
);

#if (ELMER_TIME == HARMONIC)
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

scalar maxRemDiff(readScalar(physicalProperties.lookup("maxRemDiff")));

//
// Maximum allowable relative field difference in any cell comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     >1  - once
//     1 - magnetic Reynolds number is used instead
//     (0,1) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRelDiff(readScalar(physicalProperties.lookup("maxRelDiff")));

#elif (ELMER_TIME == TRANSIENT)
int writeMultiplier( readScalar(runTime.controlDict().lookup("writeMultiplier")) );
long long int writeCounter = 0;
word writeControlDict(runTime.controlDict().lookup("writeControl"));
const bool adjustableRunTime = (writeControlDict=="adjustableRunTime");

#endif
