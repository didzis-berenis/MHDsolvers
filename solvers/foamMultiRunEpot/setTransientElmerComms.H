
// Save old U field
forAll(regionNames, i)
{
    solvers.getVelocity(regionNames[i],URegion[i]);
    vectorFieldToGlobal(UGlobal,URegion[i],regionNames[i]);
    U_oldRegion[i] = URegion[i];
}
// Send fields to Elmer
if (initialize_elmer) sending.initialize();
sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
sending.sendVector(UGlobal);
// Receive fields from Elmer
Info<< "\nReceiving fields from Elmer for all regions " << endl;
if (initialize_elmer) receiving.initialize();
receiving.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
receiving.recvVector(JGlobal);
receiving.recvVector(BGlobal);

//Lorentz force term
JxBGlobal = (JGlobal ^ BGlobal);
JJsigmaGlobal = (JGlobal & JGlobal)*sigmaGlobalInv;
//convert global field to regions
#include "updateRegionFieldsEM.H"
