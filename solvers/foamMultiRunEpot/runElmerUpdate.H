
// Save old U field
forAll(regionNames, i)
{
    if (solvers.isFluid(regionNames[i]))
    {
        const volVectorField& U = solvers.getFluid(regionNames[i]).U;
        solvers.vectorFieldToGlobal(
            UGlobal,
            U,
            regionNames[i]
        );
    }
    solvers.storeU(regionNames[i]);
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
if (solvers.isElectroHarmonic())
{
    receiving.recvVector(JimGlobal);
}
receiving.recvVector(BreGlobal);
if (solvers.isElectroHarmonic())
{
    receiving.recvVector(BimGlobal);
}
//convert global field to regions
solvers.setJ(JreGlobal);
solvers.setB(BreGlobal);
if (solvers.isElectroHarmonic())
{
    solvers.setJ(JimGlobal,true);
    solvers.setB(BimGlobal,true);
}
//Update region Lorentz force term
// and Joule heating
forAll(regionNames, i)
{
    solvers.electromagneticPredictor(regionNames[i]);
}
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
