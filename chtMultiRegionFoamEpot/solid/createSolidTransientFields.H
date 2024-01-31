// Initialise solid field pointer lists
PtrList<volVectorField> JSolid(solidRegions.size());
PtrList<volVectorField> BSolid(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "    Adding to JSolid\n" << endl;
    JSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "J",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    Info<< "    Adding to BSolid\n" << endl;
    BSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "B",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
}
