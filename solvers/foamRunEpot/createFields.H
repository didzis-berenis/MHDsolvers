// Initialize fields
//volVectorField JxB = regionSolver.getElectro().JxB;
//volScalarField JJsigma = regionSolver.getElectro().JJsigma;
volVectorField Jre = regionSolver.getElectro().J();
volVectorField Jim = regionSolver.getElectro().J(true);
volVectorField Bre = regionSolver.getElectro().B();
volVectorField Bim = regionSolver.getElectro().B(true);
const volVectorField& U = regionSolver.getU();
const volScalarField& sigma = regionSolver.getElectro().sigma();
// Cannot use const reference for Elmer functions,
// so preparing designated field for sending velocity to Elmer.
volVectorField Usent = U;
volScalarField sigmaSent = sigma;

//Initialize field names and paths
List<Pair<word>> fieldPaths = getFieldPaths(mesh);
bool needsCleanup = false;
std::set<word> fieldsToKeep = {
    "U", "p", "JxB",
    "k", "nut", "omega", "epsilon",
    "sigma", "JJsigma",
    "p_rgh", "rho",
    "PotEre","PotEim",
    "alpha.liquidmetal", "yPlus",
    "viscosity"
// Uncomment for debugging
//"sigmaInv","deltaUxB","deltaUxBre","deltaUxBim",
//"deltaJ","deltaJre","deltaJim","J","Jre","Jim","PotE",
};
// See also defaultWriteOption_ in ElectromagneticModels/electromagneticModels/electromagneticModel.H};
if (regionSolver.isFluid() || regionSolver.isSolid())//fluid == hydro + thermo
{
    //fieldsToKeep.insert("JJsigma");
    fieldsToKeep.insert("T");
}
//Viscosity field for output
volScalarField nuGlobal
(
    IOobject
    (
        "viscosity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, 0)
);

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

const dimensionedScalar absoluteZero(dimTemperature,-273.15);

