// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());
PtrList<volVectorField> JxBSolid(solidRegions.size());
PtrList<volScalarField> JJsigmaSolid(solidRegions.size());
PtrList<volVectorField> JSolid(solidRegions.size());
PtrList<volVectorField> BSolid(solidRegions.size());
PtrList<volVectorField> USolid(solidRegions.size());
PtrList<dimensionedScalar> sigmaSolid(solidRegions.size());

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

    Info<< "    Adding to JxBSolid\n" << endl;
    JxBSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "JxB",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
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

    Info<< "    Adding to USolid\n" << endl;
    USolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i]
        )
    );

    IOdictionary physicalProperties
        (
            IOobject
            (
                "physicalProperties",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    //Electrical conductivity
    sigmaSolid.set
    (
        i,
        new dimensionedScalar
        (
            "sigma",
            dimensionSet(-1, -3, 3, 0, 0, 2, 0),
            physicalProperties
        )
    );
}
