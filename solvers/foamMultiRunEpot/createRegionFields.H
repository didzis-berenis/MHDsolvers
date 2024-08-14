// Initialise field pointer lists
PtrList<volVectorField> URegion(regionNames.size());
PtrList<volVectorField> U_oldRegion(regionNames.size());
PtrList<volVectorField> JxBRegion(regionNames.size());
PtrList<volScalarField> JJsigmaRegion(regionNames.size());
PtrList<volScalarField> pRegion(regionNames.size());
PtrList<volScalarField> TRegion(regionNames.size());
PtrList<volScalarField> GeometryIdsRegion(regionNames.size());
PtrList<volScalarField> alpha1Region(regionNames.size());
List<bool> solidificationEnabled(regionNames.size(), false);
List<IOobject::writeOption> alphaWriteOption(regionNames.size(), IOobject::NO_WRITE);
PtrList<dimensionedScalar> sigmaRegion(regionNames.size());
List<double> magneticReynoldsRegion(regionNames.size());
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

    //Solidification configurations
    List<Pair<word>> regionPath = getFieldPaths(solvers.mesh(regionNames[i]));
    forAll(regionPath, j)
    {
        if (regionPath[j].first() == solverSolidificationName)
        {
            solidificationEnabled[i] = true;
            alphaWriteOption[i] = IOobject::AUTO_WRITE;
            solidificationGlobalEnabled = true;
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
    sigmaRegion.set
    (
        i,
        new dimensionedScalar
        (
            "sigma",
            dimensionSet(-1, -3, 3, 0, 0, 2, 0),
            physicalProperties
        )
    );


    #if (ELMER_TIME == HARMONIC_TIME)
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
            4*3.14159*(std::pow(10,-7))*sigmaRegion[i].value()*characteristicSize.value();
        Info<< "Rem0 = " << magneticReynoldsRegion[i] << endl;
    #endif
}
