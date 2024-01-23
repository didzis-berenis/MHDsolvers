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
    vectorMultiRegionToGlobal(UGlobal,UFluid,USolid);
    Usent = UGlobal;
}
// Send fields to Elmer
if (initialize_elmer) sending.initialize();
sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
sending.sendVector(UFluid[0]);
// Receive fields from Elmer
Info<< "\nReceiving fields from Elmer for all regions " << endl;
if (initialize_elmer) receiving.initialize();
receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
receiving.recvVector(JreGlobal);
receiving.recvVector(JimGlobal);
receiving.recvVector(BreGlobal);
receiving.recvVector(BimGlobal);

//convert global field to regions
vectorGlobalToMultiRegion(JreGlobal,JreFluid,JreSolid);
vectorGlobalToMultiRegion(BimGlobal,BimFluid,BimSolid);
forAll(fluidRegions, i)
{
    //Lorentz force term initialization
    JxBFluid[i] =  0.5*((JreFluid[i] ^ BreFluid[i])+(JimFluid[i] ^ BimFluid[i]));
    JJsigmaFluid[i] =  0.5*((JreFluid[i] & JreFluid[i])+(JimFluid[i] & JimFluid[i]))/sigmaFluid[i];
}
forAll(solidRegions, i)
{
    //Lorentz force term initialization
    JxBSolid[i] =  0.5*((JreSolid[i] ^ BreSolid[i])+(JimSolid[i] ^ BimSolid[i]));
    JJsigmaSolid[i] =  0.5*((JreSolid[i] & JreSolid[i])+(JimSolid[i] & JimSolid[i]))/sigmaSolid[i];
}
