
forAll(regionNames, i)
{
// Populate field pointer lists
    Info<< "*** Reading mesh thermophysical properties for region "
        << regionNames[i] << nl << endl;

    Info<< "    Adding to URegion\n" << endl;
    URegion.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                solvers.mesh(regionNames[i]),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solvers.mesh(regionNames[i]),
            dimensionedVector(dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::vector(0,0,0))
        )
    );
    if (solvers.isFluid(regionNames[i]))
    {
        URegion[i] = solvers.getFluid(regionNames[i]).U;
    }

    Info<< "    Adding to U_oldRegion\n" << endl;
    U_oldRegion.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U_old",
                runTime.timeName(),
                solvers.mesh(regionNames[i]),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            URegion[i]
        )
    );
/*
    Info<< "    Adding to JxBRegion\n" << endl;
    JxBRegion.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "JxB",
                runTime.timeName(),
                solvers.mesh(regionNames[i]),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(regionNames[i]),
            //not reading because solved in global mesh
            dimensionedVector(dimensionSet(1, -2, -2, 0, 0, 0, 0), Foam::vector(0,0,0))
        )
    );

    Info<< "    Adding to JJsigmaRegion\n" << endl;
    JJsigmaRegion.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "JJsigma",
                runTime.timeName(),
                solvers.mesh(regionNames[i]),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(regionNames[i]),
            //not reading because solved in global mesh
            dimensionedScalar(dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
        )
    );

*/
    Info<< "    Adding to pRegion\n" << endl;
    pRegion.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.constant(),
                solvers.mesh(regionNames[i]),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(regionNames[i]),
            dimensionedScalar(dimPressure, 0)
        )
    );

    Info<< "    Adding to TRegion\n" << endl;
    TRegion.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.constant(),
                solvers.mesh(regionNames[i]),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(regionNames[i]),
            dimensionedScalar(dimTemperature, 0)
        )
    );

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

    IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                solvers.mesh(regionNames[i]),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    //Electrical conductivity
    /*sigmaRegionConst.set
    (
        i,
        new dimensionedScalar
        (
            "sigma",
            dimensionSet(-1, -3, 3, 0, 0, 2, 0),
            physicalProperties
        )
    );*/
    dimensionedScalar sigmaRegionConst
    (
        "sigma",
        dimensionSet(-1, -3, 3, 0, 0, 2, 0),
        physicalProperties
    );

    if (solvers.isElectroHarmonic())
    {
    //#if (ELMER_TIME == HARMONIC_TIME)
        //Characteristic length
        dimensionedScalar characteristicSize
        (
            "Lchar",
            dimensionSet(0, 1, 0, 0, 0, 0, 0),
            0
        );
        
        if (solvers.isFluid(regionNames[i]))
        {
            characteristicSize = dimensionedScalar
            (
                "Lchar",
                dimensionSet(0, 1, 0, 0, 0, 0, 0),
                physicalProperties
            );
        }

        magneticReynoldsRegion[i] = 
            4*3.14159*(std::pow(10,-7))*sigmaRegionConst.value()*characteristicSize.value();
        Info<< "Rem0 = " << magneticReynoldsRegion[i] << endl;
    }//#endif
}
