
// Save old U field
forAll(fluidRegions, i)
{
    U_oldFluid[i] = UFluid[i];
}
//Use only fluid mesh if single fluid region
if (fluidRegions.size() == 1)
{
    Info<< "\nSending fields to Elmer for a single region " << endl;
    Usent = UFluid[0];
}
//Use global mesh if multiple fluid regions
else
{ 
    Info<< "\nSending fields to Elmer for all regions " << endl;
    //convert global field to regions
    forAll(fluidRegions, i)
    {
        vectorFluidFieldToGlobal(UGlobal,Usent,i);
    }
    Usent = UGlobal;
}
// Send fields to Elmer
if (initialize_elmer) sending.initialize();
sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
sending.sendVector(UFluid[0]);
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
vectorGlobalToMultiRegion(JxBGlobal,JxBFluid,JxBSolid);
scalarGlobalToMultiRegion(JJsigmaGlobal,JJsigmaFluid,JJsigmaSolid);
