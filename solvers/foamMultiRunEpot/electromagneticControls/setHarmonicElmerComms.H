
// Save old U field
/*forAll(regionNames, i)
{
    if (solvers.isFluid(regionNames[i]))
    {
        URegion[i] = solvers.getVelocity(regionNames[i]);
    }
    vectorFieldToGlobal(UGlobal,URegion[i],regionNames[i]);
    U_oldRegion[i] = URegion[i];
}*/
UGlobal = solvers.getGlobalVelocity();
//vectorMultiRegionToGlobal(U_oldGlobal,U_oldRegion);
// Send fields to Elmer
if (initialize_elmer) sending.initialize();
sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
sending.sendVector(UGlobal);
// Receive fields from Elmer
Info<< "\nReceiving fields from Elmer for all regions " << endl;
if (initialize_elmer) receiving.initialize();
receiving.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
receiving.recvVector(JreGlobal);
receiving.recvVector(JimGlobal);
receiving.recvVector(BreGlobal);
receiving.recvVector(BimGlobal);

solvers.updateRegionJ(JreGlobal);
solvers.updateRegionB(BreGlobal);
solvers.updateRegionJ(JimGlobal,true);
solvers.updateRegionB(BimGlobal,true);

//Lorentz force term
JxBGlobal =  0.5*((JreGlobal ^ BreGlobal)+(JimGlobal^ BimGlobal));
JJsigmaGlobal =  0.5*((JreGlobal & JreGlobal)+(JimGlobal & JimGlobal))*sigmaGlobalInv;
//convert global field to regions
//#include "updateRegionFieldsEM.H"
