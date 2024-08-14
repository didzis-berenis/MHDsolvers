//
// EOF-Library supports using single mesh of communication either for sending or receiving
//
Info<< "*** Preparing Elmer communications for sending" << nl << endl;
Elmer<fvMesh> sending(
        meshGlobal, //mesh
        1, // 1=send, -1=receive
        0, // 1=initialize, 0=w/o init
        1 // 1=multiregion, 0=exports O2E files
    );

Info<< "*** Preparing Elmer communications for receiving" << nl << endl;
Elmer<fvMesh> receiving(
        meshGlobal, //mesh
        -1, // 1=send, -1=receive
        0, // 1=initialize, 0=w/o init
        1 // 1=multiregion, 0=exports O2E files
    );

// Get control parameters
//#if (ELMER_TIME == HARMONIC_TIME)
if (solvers.isElectroHarmonic())
{
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
}
else
//#elif (ELMER_TIME == TRANSIENT_TIME)
{
    int writeMultiplier( readScalar(runTime.controlDict().lookup("writeMultiplier")) );
    long long int writeCounter = 0;
    word writeControlDict(runTime.controlDict().lookup("writeControl"));
    const bool adjustableRunTime = (writeControlDict=="adjustableRunTime");

}
//#endif