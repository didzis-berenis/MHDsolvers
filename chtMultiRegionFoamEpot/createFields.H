#include "createFluidFields.H"
#include "createSolidFields.H"
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
