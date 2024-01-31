// Initialise solid field pointer lists
PtrList<volVectorField> JreSolid(solidRegions.size());
PtrList<volVectorField> BreSolid(solidRegions.size());
PtrList<volVectorField> JimSolid(solidRegions.size());
PtrList<volVectorField> BimSolid(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "    Adding to JreSolid\n" << endl;
    JreSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Jre",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    Info<< "    Adding to JimSolid\n" << endl;
    JimSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Jim",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    Info<< "    Adding to BreSolid\n" << endl;
    BreSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Bre",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    Info<< "    Adding to BimSolid\n" << endl;
    BimSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Bim",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
}
