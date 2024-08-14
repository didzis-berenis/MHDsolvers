// Initialise field pointer lists
PtrList<volVectorField> URegion(regionNames.size());
PtrList<volVectorField> U_oldRegion(regionNames.size());
//PtrList<volVectorField> JxBRegion(regionNames.size());
//PtrList<volScalarField> JJsigmaRegion(regionNames.size());
PtrList<volScalarField> pRegion(regionNames.size());
PtrList<volScalarField> TRegion(regionNames.size());
PtrList<volScalarField> GeometryIdsRegion(regionNames.size());
PtrList<volScalarField> alpha1Region(regionNames.size());
//PtrList<dimensionedScalar> sigmaRegionConst(regionNames.size());

List<double> magneticReynoldsRegion(regionNames.size());

//Solidification configurations
word solverSolidificationName = "solidification:alpha1";
word userSolidificationName = "alpha1";
bool solidificationGlobalEnabled = false;
List<bool> solidificationEnabled(regionNames.size(), false);
List<IOobject::writeOption> alphaWriteOption(regionNames.size(), IOobject::NO_WRITE);

//Initialize field names and paths
List<List<Pair<word>>> regionPaths(regionNames.size());

std::set<word> fieldsToKeep = {"U", "p", "JxB", "JJsigma", "T"};//, "GeometryIds"

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
/*
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
*/
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

Info<< "Initializing field JreGlobal\n" << endl;
volVectorField JreGlobal
(
    IOobject
    (
        "Jre",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field JimGlobal\n" << endl;
volVectorField JimGlobal
(
    IOobject
    (
        "Jim",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field BreGlobal\n" << endl;
volVectorField BreGlobal
(
    IOobject
    (
        "Bre",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
Info<< "Initializing field BimGlobal\n" << endl;
volVectorField BimGlobal
(
    IOobject
    (
        "Bim",
        runTime.timeName(),
        meshGlobal,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshGlobal
);
/*
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
*/
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
forAll(regionNames, i)
{
    const volScalarField& sigmaInvRegion = solvers.getElectro(regionNames[i]).sigmaInv();
    solvers.scalarFieldToGlobal(
        sigmaGlobalInv,
        sigmaInvRegion,
        regionNames[i]
    );
}
/*
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
*/
/*
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
*/
//List<Pair<word>> fieldPaths = getFieldPaths(meshGlobal);

//auto alphaWriteOptionGlobal = IOobject::NO_WRITE;
