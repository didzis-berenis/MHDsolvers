// Initialise field pointer lists
PtrList<volVectorField> U_oldRegion(regionSolverNames.size());
PtrList<volVectorField> JxBRegion(regionSolverNames.size());
PtrList<volScalarField> JJsigmaRegion(regionSolverNames.size());
PtrList<volScalarField> GeometryIdsRegion(regionSolverNames.size());
PtrList<dimensionedScalar> sigmaRegion(regionSolverNames.size());
List<double> magneticReynoldsRegion(regionSolverNames.size());
forAll(regionSolverNames, i)
{
// Populate field pointer lists
    Info<< "*** Reading mesh thermophysical properties for region "
        << regionSolverNames[i].first() << nl << endl;
        
    const word& solverName = regionSolverNames[i].second();

    Info<< "    Adding to U_oldRegion\n" << endl;
    if (solverName == "conductingFluid")
    {
        U_oldRegion.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U_old",
                    runTime.timeName(),
                    solvers.mesh(i),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solvers.fluids(i)->getVelocity()
            )
        );
    }
    else
    {
        U_oldRegion.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U_old",
                    runTime.timeName(),
                    solvers.mesh(i),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solvers.mesh(i),
                dimensionedVector(dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::vector(0,0,0))
            )
        );

    }

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
                solvers.mesh(i),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(i),
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
                solvers.mesh(i),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(i),
            //not reading because solved in global mesh
            dimensionedScalar(dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
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
                solvers.mesh(i),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solvers.mesh(i)
        )
    );

    IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                solvers.mesh(i),
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
        
        if (solverName == "conductingFluid")
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
