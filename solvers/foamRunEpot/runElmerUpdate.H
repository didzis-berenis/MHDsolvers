
// Save old U field
regionSolver.storeU();

// Save old phase-fraction field
regionSolver.storeAlpha1();

// Send fields to Elmer
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
Info<< "\nReceiving fields from Elmer " << endl;
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
//update solver fields
regionSolver.setJ(Jre);
regionSolver.setB(Bre);
if (regionSolver.isElectroHarmonic())
{
    regionSolver.setJ(Jim,true);
    regionSolver.setB(Bim,true);
}
//Update mixture conductivity if VoF is used
regionSolver.updateMixtureConductivity();
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
