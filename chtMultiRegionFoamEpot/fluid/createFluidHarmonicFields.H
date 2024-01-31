// Initialise fluid field pointer lists
PtrList<volVectorField> JreFluid(fluidRegions.size());
PtrList<volVectorField> BreFluid(fluidRegions.size());
PtrList<volVectorField> JimFluid(fluidRegions.size());
PtrList<volVectorField> BimFluid(fluidRegions.size());
PtrList<volScalarField> PotEreFluid(fluidRegions.size());
PtrList<volScalarField> PotEimFluid(fluidRegions.size());
List<double> magneticReynolds(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "    Adding to JreFluid\n" << endl;
    JreFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Jre",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to JimFluid\n" << endl;
    JimFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Jim",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to BreFluid\n" << endl;
    BreFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Bre",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to BimFluid\n" << endl;
    BimFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Bim",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    //Electric potential field
    Info<< "    Adding to PotEreFluid\n" << endl;
    PotEreFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "PotEre",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to PotEimFluid\n" << endl;
    PotEimFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "PotEim",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    fluidRegions[i].schemes().setFluxRequired(PotEreFluid[i].name());

    fluidRegions[i].schemes().setFluxRequired(PotEimFluid[i].name());

    IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    //Characteristic length
    dimensionedScalar characteristicSize 
        (
            "Lchar",
            dimensionSet(0, 1, 0, 0, 0, 0, 0),
            physicalProperties
        );

	magneticReynolds[i] = 
        4*3.14159*(std::pow(10,-7))*sigmaFluid[i].value()*characteristicSize.value();
    Info<< "Rem0 = " << magneticReynolds[i] << endl;
}
