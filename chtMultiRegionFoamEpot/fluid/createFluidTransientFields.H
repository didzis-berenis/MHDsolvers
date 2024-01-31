// Initialise fluid field pointer lists
PtrList<volVectorField> JFluid(fluidRegions.size());
PtrList<volVectorField> BFluid(fluidRegions.size());
PtrList<volScalarField> PotEFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "    Adding to JFluid\n" << endl;
    JFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "J",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to BFluid\n" << endl;
    BFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "B",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    //Electric potential field
    Info<< "    Adding to PotEFluid\n" << endl;
    PotEFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "PotE",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    fluidRegions[i].schemes().setFluxRequired(PotEFluid[i].name());
}
