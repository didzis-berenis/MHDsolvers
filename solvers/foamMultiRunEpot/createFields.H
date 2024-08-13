//Solidification configurations
word solverSolidificationName = "solidification:alpha1";
word userSolidificationName = "alpha1";
bool solidificationGlobalEnabled = false;

#include "createRegionFields.H"
volVectorField UGlobal
(
    IOobject
    (
        "U",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::vector(0,0,0))
);
volVectorField U_oldGlobal
(
    IOobject
    (
        "U_old",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    UGlobal
);
volScalarField pGlobal
(
    IOobject
    (
        "p",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimPressure, 0)
);
volScalarField TGlobal
(
    IOobject
    (
        "T",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimTemperature, 0)
);
volScalarField JJsigmaGlobal
(
    IOobject
    (
        "JJsigma",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
);

volVectorField JxBGlobal
(
    IOobject
    (
        "JxB",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, -2, -2, 0, 0, 0, 0), Foam::vector(0,0,0))
);
volScalarField GeometryIdsGlobal
(
    IOobject
    (
        "GeometryIds",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
);
volScalarField sigmaGlobal
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(-1, -3, 3, 0, 0, 2, 0), 1)
);
volScalarField sigmaGlobalInv
(
    IOobject
    (
        "sigmaInv",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(1, 3, -3, 0, 0, -2, 0), 1)
);
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(1, 0, -1, 0, 0, 0, 0), 0)
);
scalarMultiRegionToGlobal(GeometryIdsGlobal,GeometryIdsRegion);

forAll(sigmaGlobal, cellI)
{
    scalar sigma = sigmaRegion[cellToRegion[cellI]].value();
    sigmaGlobal[cellI] = sigma == 0 ? ROOTVSMALL : sigma;
    sigmaGlobalInv[cellI] = sigma == 0 ? 0 : 1/sigma;
}

#if (ELMER_TIME == HARMONIC_TIME)
    #include "createHarmonicFields.H"
#elif (ELMER_TIME == TRANSIENT_TIME)
    #include "createTransientFields.H"
#endif

//Initialize field names and paths
List<word> regionPaths(regionNames.size());

forAll(regionNames, i)
{
    regionPaths[i] = getRegionPath(solvers.mesh(regionNames[i]));
}
//List<Pair<word>> fieldPaths = getFieldPaths(meshGlobal);
std::set<word> fieldsToKeep = {"U", "p", "JxB", "JJsigma", "T", "GeometryIds"};

auto alphaWriteOptionGlobal = IOobject::NO_WRITE;
if (solidificationGlobalEnabled)
{
    fieldsToKeep.insert(userSolidificationName);
    alphaWriteOptionGlobal = IOobject::AUTO_WRITE;
}
volScalarField alpha1Global
(
    IOobject
    (
        userSolidificationName,
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        alphaWriteOptionGlobal
    ),
    meshGlobal,
    dimensionedScalar(dimless, 0)
);
if (solidificationGlobalEnabled)
{
	scalarMultiRegionToGlobal(alpha1Global,alpha1Region);
}

//Update field paths to include alpha1
//fieldPaths = getFieldPaths(meshGlobal);
/*
forAll(fieldPaths, i)
{
    word fieldName = fieldPaths[i].first();
    keepField[fieldName] = fieldsToKeep.find(fieldName) != fieldsToKeep.end();
}*/

