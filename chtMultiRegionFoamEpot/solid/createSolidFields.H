// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());
PtrList<volScalarField> JJsigmaSolid(fluidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoSolid\n" << endl;
    thermoSolid.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );

    Info<< "    Adding to JJsigmaSolid\n" << endl;
    JJsigmaSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "JJsigma",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
}
