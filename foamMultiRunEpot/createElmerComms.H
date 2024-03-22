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

#if (ELMER_TIME == HARMONIC_TIME)
    #include "createHarmonicElmerComms.H"
#elif (ELMER_TIME == TRANSIENT_TIME)
    #include "createTransientElmerComms.H"
#endif