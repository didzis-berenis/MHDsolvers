// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());
PtrList<volVectorField> JxBSolid(solidRegions.size());
PtrList<volScalarField> JJsigmaSolid(solidRegions.size());
PtrList<volScalarField> GeometryIdsSolid(solidRegions.size());
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
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            //not reading because solved in global mesh
            dimensionedVector(dimensionSet(1, -2, -2, 0, 0, 0, 0), Foam::vector(0,0,0))
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
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            //not reading because solved in global mesh
            dimensionedScalar(dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
        )
    );

    Info<< "    Adding to GeometryIdsSolid\n" << endl;
    GeometryIdsSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "GeometryIds",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
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
