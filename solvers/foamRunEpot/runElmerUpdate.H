
// Save old U field
regionSolver.storeU();

// Save old phase-fraction field
regionSolver.storeAlpha1();

for (int iter = 0; iter < nElmerCorrectors; iter++)
{
    volVectorField pJre = Jre;
    volVectorField pJim = Jim;
    volVectorField pBre = Bre;
    volVectorField pBim = Bim;
    Info << "\nElmer iteration: " << iter << endl;
    // Send fields to Elmer
    Info << "Sending fields to Elmer" << endl;
    if (initialize_elmer) sending.initialize();
    sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
    //cannot use const reference, so assign to new field 
    Usent = U;
    sigmaSent = sigma;
    sending.sendVector(Usent);
    if (regionSolver.isIncompressibleConductingVoF())
    {
        sending.sendScalar(sigmaSent);
    }
    // Receive fields from Elmer
    Info << "Receiving fields from Elmer" << endl;
    if (initialize_elmer) receiving.initialize();
    receiving.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(Jre);
    if (regionSolver.isElectroHarmonic())
    {
        receiving.recvVector(Jim);
    }
    receiving.recvVector(Bre);
    if (regionSolver.isElectroHarmonic())
    {
        receiving.recvVector(Bim);
    }
    scalar eJre = gMax(mag(pJre - Jre)())/(gMax(mag(pJre)())+SMALL);
    scalar eJim = gMax(mag(pJim - Jim)())/(gMax(mag(pJim)())+SMALL);
    scalar eBre = gMax(mag(pBre - Bre)())/(gMax(mag(pBre)())+SMALL);
    scalar eBim = gMax(mag(pBim - Bim)())/(gMax(mag(pBim)())+SMALL);
    scalar maxError = max(max(eJre,eJim),max(eBre,eBim));
    bool quitLoop = maxError < ROOTSMALL;
    // Initialization is needed only on the first update.
    if (initialize_elmer)
        initialize_elmer = false;
    //Quit the loop if bad status
    quitLoop |= elmer_status != 1;
    if (quitLoop) break;
}
//update solver fields
regionSolver.setJ(Jre);
regionSolver.setB(Bre);
if (regionSolver.isElectroHarmonic())
{
    regionSolver.setJ(Jim,true);
    regionSolver.setB(Bim,true);
}
//Update solver Lorentz force term
// and Joule heating
regionSolver.electromagneticPredictor();
			
// Log the current simulation time
if (Pstream::master())
{
    std::ofstream elmerTimes(elmerTimesFileName, std::ios::app);
    if (elmerTimes.is_open())
    {
        elmerTimes << runTime.timeName() << std::endl;
        elmerTimes.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << elmerTimesFileName << " for writing!\n" << abort(FatalError);
}
