/*
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

//Info<< "Initializing field Jre\n" << endl;
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

//Info<< "Initializing field Jim\n" << endl;
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

//Info<< "Initializing field Bre\n" << endl;
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

//Info<< "Initializing field Bim\n" << endl;
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
*/
volVectorField JxB = regionSolver.getElectro().JxB;
volScalarField JJsigma = regionSolver.getElectro().JJsigma;
volVectorField Jre = regionSolver.getElectro().J();
volVectorField Jim = regionSolver.getElectro().J(true);
volVectorField Bre = regionSolver.getElectro().B();
volVectorField Bim = regionSolver.getElectro().B(true);
const volScalarField& sigmaInv = regionSolver.getElectro().sigmaInv();
const volVectorField& U = regionSolver.getVelocity();
// Cannot use const reference for Elmer functions,
// so preparing designated field for sending velocity to Elmer.
volVectorField Usent = U;

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
