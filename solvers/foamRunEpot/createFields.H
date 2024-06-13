
volScalarField& p = regionSolver.getPressure();

Info<< "Initializing field U\n" << endl;
volVectorField& U = regionSolver.getVelocity();

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
    0
);

if (regionSolver.isFluid())
{
    sigma = dimensionedScalar
    (
        "sigma",
        dimensionSet(-1, -3, 3, 0, 0, 2, 0),
        physicalProperties
    );
}

Info<< "Initializing field JxB\n" << endl;
volVectorField JxB
(
    IOobject
    (
        "JxB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(1, -2, -2, 0, 0, 0, 0), Foam::vector(0,0,0))
);

Info<< "Initializing field JJsigma\n" << endl;
volScalarField JJsigma
(
    IOobject
    (
        "JJsigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
);

#if (ELMER_TIME == HARMONIC_TIME)
    #include "createHarmonicFields.H"
#elif (ELMER_TIME == TRANSIENT_TIME)
    #include "createTransientFields.H"
#endif

//Initialize field names and paths
List<Pair<word>> fieldPaths = getFieldPaths(mesh);
std::set<word> fieldsToKeep = {"U", "p", "JxB", "JJsigma", "T"};

//Solidification configurations
word solverSolidificationName = "solidification:alpha1";
word userSolidificationName = "alpha1";
bool solidificationEnabled = false;
auto alphaWriteOption = IOobject::NO_WRITE;
forAll(fieldPaths, i)
{
    if (fieldPaths[i].first() == solverSolidificationName)
    {
        solidificationEnabled = true;
        fieldsToKeep.insert(userSolidificationName);
        alphaWriteOption = IOobject::AUTO_WRITE;
    }
}
volScalarField alpha1
(
    IOobject
    (
        userSolidificationName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        alphaWriteOption
    ),
    mesh,
    dimensionedScalar(dimless, 0)
);
if (solidificationEnabled)
{
	alpha1 = mesh().lookupObject<volScalarField>(solverSolidificationName);
}

//Update field paths to include alpha1
fieldPaths = getFieldPaths(mesh);

forAll(fieldPaths, i)
{
    word fieldName = fieldPaths[i].first();
    keepField[fieldName] = fieldsToKeep.find(fieldName) != fieldsToKeep.end();
}
