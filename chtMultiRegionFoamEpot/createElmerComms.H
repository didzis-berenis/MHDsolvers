// Initialise Elmer sending pointer lists for fluid regions
PtrList<Elmer<fvMesh>> sendingFluid(fluidRegions.size());
PtrList<Elmer<fvMesh>> receivingFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Preparing Elmer communications for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to sendingFluid\n" << endl;
    sendingFluid.set
    (
        i,
        new Elmer<fvMesh>
        (
            fluidRegions[i], //mesh
            1, // 1=send, -1=receive
            0, // 1=initialize, 0=w/o init
            1 // 1=multiregion, 0=exports O2E files
        )
    );

    Info<< "    Adding to receivingFluid\n" << endl;
    receivingFluid.set
    (
        i,
        new Elmer<fvMesh>
        (
            fluidRegions[i], //mesh
            -1, // 1=send, -1=receive
            0, // 1=initialize, 0=w/o init
            1 // 1=multiregion, 0=exports O2E files
        )
    );

}

// Initialise Elmer sending pointer lists for solid regions
PtrList<Elmer<fvMesh>> receivingSolid(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Preparing Elmer communications for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to receivingSolid\n" << endl;
    receivingSolid.set
    (
        i,
        new Elmer<fvMesh>
        (
            solidRegions[i], //mesh
            -1, // 1=send, -1=receive
            0, // 1=initialize, 0=w/o init
            1 // 1=multiregion, 0=exports O2E files
        )
    );

}

// Maximum allowable magnetic Reynolds number difference comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     (0,inf) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRemDiff(readScalar(runTime.controlDict().lookup("maxRemDiff")));

//
// Maximum allowable relative field difference in any cell comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     >1  - once
//     1 - magnetic Reynolds number is used instead
//     (0,1) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRelDiff(readScalar(runTime.controlDict().lookup("maxRelDiff")));
