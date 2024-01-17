// Send fields to Elmer
forAll(fluidRegions, i)
{
    U_oldFluid[i] = UFluid[i];
    Info<< "\nSending fields to Elmer for fluid region "
        << fluidRegions[i].name() << endl;
    Elmer<fvMesh>& sending = sendingFluid[i];
    if (initialize_elmer) sending.initialize();
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    sending.sendVector(UFluid[i]);
}

Info<< "\nReceiving fields from Elmer for all regions " << endl;
if (initialize_elmer) receivingGlobal.initialize();
receivingGlobal.sendStatus(1); // 1=ok, 0=lastIter, -1=error
receivingGlobal.recvVector(JGlobal);
receivingGlobal.recvVector(BGlobal);
/*
//convert global field to regions
globalToLocal(JGlobal,JFluid,JSolid);
globalToLocal(BGlobal,BFluid,BSolid);
forAll(fluidRegions, i)
{
    //Lorentz force term initialization
    JxBFluid[i] =  (JFluid[i] ^ BFluid[i]);
    dimensionedScalar& sigma = sigmaFluid[i];
    JJsigmaFluid[i] =  (JFluid[i] & JFluid[i])/sigma;
}
forAll(solidRegions, i)
{
    //Lorentz force term initialization
    JxBSolid[i] =  (JSolid[i] ^ BSolid[i]);
    dimensionedScalar& sigma = sigmaSolid[i];
    JJsigmaSolid[i] =  (JSolid[i] & JSolid[i])/sigma;
}
*/

