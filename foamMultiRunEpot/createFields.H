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
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(1, 3, -3, 0, 0, -2, 0), 1)
);
//0=Write fields separately for each region
//1=Write fields in a merged(global) mesh
const bool writeMerged( readScalar(runTime.controlDict().lookup("writeMerged")) );
scalarMultiRegionToGlobal(GeometryIdsGlobal,GeometryIdsRegion);
List<scalar> conductivity(regionSolverNames.size());

forAll(regionSolverNames, i)
{
    const word& regionName = regionSolverNames[i].first();
    //const word& solverName = regionSolverNames[i].second();
    conductivity[regionID[regionName]] = sigmaRegion[i].value();
}
forAll(sigmaGlobal, cellI)
{
    scalar sigma = conductivity[cellToRegion[cellI]];
    sigmaGlobal[cellI] = sigma == 0 ? ROOTVSMALL : sigma;
    sigmaGlobalInv[cellI] = sigma == 0 ? 0 : 1/sigma;
}

#if (ELMER_TIME == HARMONIC_TIME)
    #include "createHarmonicFields.H"
#elif (ELMER_TIME == TRANSIENT_TIME)
    #include "createTransientFields.H"
#endif
