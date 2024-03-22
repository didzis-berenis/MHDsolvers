
// Save old U field
forAll(regionSolverNames, i)
{
    const word& solverName = regionSolverNames[i].second();
    if (solverName == "conductingFluid")
    {
        U_oldRegion[i] = solvers.fluids(fluidID[solverName])->getVelocity();
    } 
}
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

//Lorentz force term
JxBGlobal =  0.5*((JreGlobal ^ BreGlobal)+(JimGlobal^ BimGlobal));
JJsigmaGlobal =  0.5*((JreGlobal & JreGlobal)+(JimGlobal & JimGlobal))*sigmaGlobalInv;
//convert global field to regions
#include "updateRegionFieldsEM.H"
