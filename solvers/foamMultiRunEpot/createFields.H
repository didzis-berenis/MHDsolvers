#include "createRegionFields.H"
volVectorField UGlobal
(
    IOobject
    (
        "U",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
    //read to get boundary conditions
    //that will be set before writing global field
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
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
    //read to get boundary conditions
    //that will be set before writing global field
);
volScalarField TGlobal
(
    IOobject
    (
        "T",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
    //read to get boundary conditions
    //that will be set before writing global field
);
volScalarField JJsigmaGlobal
(
    IOobject
    (
        "JJsigma",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);

volVectorField JxBGlobal
(
    IOobject
    (
        "JxB",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
volScalarField GeometryIdsGlobal
(
    IOobject
    (
        "GeometryIds",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
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

List<word> regionPaths(regionNames.size());
forAll(regionNames, i)
{
    regionPaths[i] = getRegionPath(solvers.mesh(regionNames[i]));
}

List<Pair<word>> fieldPaths = getFieldPaths(meshGlobal);
std::set<word> fieldsToKeep = {"U", "p", "JxB", "JJsigma", "T", "GeometryIds"};
forAll(fieldPaths, i)
{
    word fieldName = fieldPaths[i].first();
    keepField[fieldName] = fieldsToKeep.find(fieldName) != fieldsToKeep.end();
}
