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

// Receive fields from Elmer
forAll(fluidRegions, i)
{
    Info<< "\nReceiving fields from Elmer for fluid region "
        << fluidRegions[i].name() << endl;
    Elmer<fvMesh>& receiving = receivingFluid[i];
    if (initialize_elmer) receiving.initialize();
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(JreFluid[i]);
    receiving.recvVector(JimFluid[i]);
    receiving.recvVector(BreFluid[i]);
    receiving.recvVector(BimFluid[i]);

    //Lorentz force term initialization
    JxBFluid[i] =  0.5*((JreFluid[i] ^ BreFluid[i]) + (JimFluid[i] ^ BimFluid[i]) );
    dimensionedScalar& sigma = sigmaFluid[i];
    JJsigmaFluid[i] =  0.5*((JreFluid[i] & JreFluid[i]) + (JimFluid[i] & JimFluid[i]) )/sigma;
}
forAll(solidRegions, i)
{
    Info<< "\nReceiving fields from Elmer for solid region "
        << solidRegions[i].name() << endl;
    Elmer<fvMesh>& receiving = receivingSolid[i];
    if (initialize_elmer) receiving.initialize();
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvScalar(JJsigmaSolid[i]);
}
