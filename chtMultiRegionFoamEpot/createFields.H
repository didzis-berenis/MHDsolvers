#include "createFluidFields.H"
#include "createSolidFields.H"
//0=Write fields separately for each region
//1=Write fields in a merged(global) mesh
const bool writeMerged( readScalar(runTime.controlDict().lookup("writeMerged")) );
volVectorField UGlobal
(
    IOobject
    (
        "U",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::vector(0,0,0))
);
volScalarField TGlobal
(
    IOobject
    (
        "T",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(0, 0, 0, 1, 0, 0, 0), 0)
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
scalarMultiRegionToGlobal(GeometryIdsGlobal,GeometryIdsFluid,GeometryIdsSolid);
volScalarField pGlobal
(
    IOobject
    (
        "p",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(1, -1, -2, 0, 0, 0, 0), 0)
);
volScalarField JJsigmaGlobal
(
    IOobject
    (
        "JJsigma",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal,
    dimensionedScalar(dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
);
//Use only fluid mesh if single fluid region
//Use global mesh if multiple fluid regions
volVectorField Usent = fluidRegions.size() == 1 ? UFluid[0] : UGlobal;

#if (ELMER_TIME == HARMONIC_TIME)
    #include "createFluidHarmonicFields.H"
    #include "createSolidHarmonicFields.H"
    #include "createHarmonicFields.H"
#elif (ELMER_TIME == TRANSIENT_TIME)
    #include "createFluidTransientFields.H"
    #include "createSolidTransientFields.H"
    #include "createTransientFields.H"
#endif
