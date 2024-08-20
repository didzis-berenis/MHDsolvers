
Info<< "Initializing field UGlobal\n" << endl;
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
Info<< "Initializing field JJsigmaGlobal\n" << endl;
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
Info<< "Initializing field JxBGlobal\n" << endl;
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

volVectorField JreGlobal
(
    IOobject
    (
        "Jre",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(0, -2, 0, 0, 0, 1, 0), Foam::vector(0,0,0))
);
volVectorField JimGlobal
(
    IOobject
    (
        "Jim",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(0, -2, 0, 0, 0, 1, 0), Foam::vector(0,0,0))
);
volVectorField BreGlobal
(
    IOobject
    (
        "Bre",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, 0, -2, 0, 0, -1, 0), Foam::vector(0,0,0))
);
volVectorField BimGlobal
(
    IOobject
    (
        "Bim",
        runTime.timeName(),
        meshGlobal,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    meshGlobal,
    dimensionedVector(dimensionSet(1, 0, -2, 0, 0, -1, 0), Foam::vector(0,0,0))
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

//Region fields
PtrList<volScalarField> GeometryIdsRegion(regionNames.size());
PtrList<volScalarField> alpha1Region(regionNames.size());
//Solidification configurations
word solverSolidificationName = "solidification:alpha1";
word userSolidificationName = "alpha1";
//bool solidificationGlobalEnabled = false;
List<bool> solidificationEnabled(regionNames.size(), false);
List<IOobject::writeOption> alphaWriteOption(regionNames.size(), IOobject::NO_WRITE);

//Initialize field names and paths
List<List<Pair<word>>> regionPaths(regionNames.size());

std::set<word> fieldsToKeep = {"U", "p", "JxB", "JJsigma", "T", "GeometryIds"};//

forAll(regionNames, i)
{
// Populate field pointer lists
    Info<< "*** Reading mesh thermophysical properties for region "
        << regionNames[i] << nl << endl;

    Info<< "    Adding to GeometryIdsRegion\n" << endl;
    GeometryIdsRegion.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "GeometryIds",
                runTime.timeName(),
                solvers.mesh(regionNames[i]),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(regionNames[i])
        )
    );

    regionPaths[i] = getFieldPaths(solvers.mesh(regionNames[i]));
    forAll(regionPaths[i], j)
    {
        word fieldName = regionPaths[i][j].first();
        keepField[Pair<word>(fieldName,regionNames[i])] = fieldsToKeep.find(fieldName) != fieldsToKeep.end();
        if (fieldName == solverSolidificationName)
        {
            keepField[Pair<word>(fieldName,regionNames[i])] = true;
            solidificationEnabled[i] = true;
            alphaWriteOption[i] = IOobject::AUTO_WRITE;
            //solidificationGlobalEnabled = true;
        }
    }

    alpha1Region.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                userSolidificationName,
                runTime.constant(),
                solvers.mesh(regionNames[i]),
                IOobject::NO_READ,
                alphaWriteOption[i]
            ),
            solvers.mesh(regionNames[i]),
            dimensionedScalar(dimless, 0)
        )
    );
    //Get alpha1 field from registry
    if (solidificationEnabled[i])
    {
        alpha1Region[i] = solvers.mesh(regionNames[i]).lookupObject<volScalarField>(solverSolidificationName);
    }
    const volScalarField& sigmaInvRegion = solvers.getElectro(regionNames[i]).sigmaInv();
    solvers.scalarFieldToGlobal(
        sigmaGlobalInv,
        sigmaInvRegion,
        regionNames[i]
    );
}
/*
scalarMultiRegionToGlobal(GeometryIdsGlobal,GeometryIdsRegion);

forAll(sigmaGlobal, cellI)
{
    scalar sigma = sigmaRegion[cellToRegion[cellI]].value();
    sigmaGlobal[cellI] = sigma == 0 ? ROOTVSMALL : sigma;
    sigmaGlobalInv[cellI] = sigma == 0 ? 0 : 1/sigma;
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
*/
