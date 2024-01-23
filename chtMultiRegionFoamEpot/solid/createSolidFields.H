// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());
PtrList<volVectorField> JxBSolid(solidRegions.size());
PtrList<volScalarField> JJsigmaSolid(solidRegions.size());
PtrList<volVectorField> JreSolid(solidRegions.size());
PtrList<volVectorField> BreSolid(solidRegions.size());
PtrList<volVectorField> JimSolid(solidRegions.size());
PtrList<volVectorField> BimSolid(solidRegions.size());
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
