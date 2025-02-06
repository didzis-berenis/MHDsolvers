
// Clean global current fields before applying sources
JreGlobal = Zero;
if (solvers.isElectroHarmonic())
{
    JimGlobal = Zero;
}
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
    if (solvers.isElectric(regionNames[i]) && solvers.getElectro(regionNames[i]).isSource())
    {
        // For these regions current density is calculated on OpenFOAM side
        // and incorporated in Elmer as an external current source. 
        //Info<< "\nAdding to global J: " << regionNames[i] << endl;
        hasElectricSources = true;
        const volVectorField& Jre = solvers.getElectro(regionNames[i]).J();
        solvers.vectorFieldToGlobal(
            JreGlobal,
            Jre,
            regionNames[i]
        );
        if (solvers.isElectroHarmonic())
        {
            const volVectorField& Jim = solvers.getElectro(regionNames[i]).J(true);
            solvers.vectorFieldToGlobal(
                JimGlobal,
                Jim,
                regionNames[i]
            );
        }
    }
}
if (transient_em)
{
    JreGlobal *= Foam::sin(omega_em*runTime.userTimeValue());
}
//JreGlobal.write();
// Send fields to Elmer
Info<< "\nSending fields to Elmer for all regions " << endl;
if (initialize_elmer) sending.initialize();
sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
sending.sendVector(UGlobal);
if (hasElectricSources)
{
    sending.sendVector(JreGlobal);
    if (solvers.isElectroHarmonic())
    {
        sending.sendVector(JimGlobal);
    }
}
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


forAll(regionNames, i)
{
    // Skip update for electric sources, since current was calculated in OpnFOAM
    if (solvers.isElectric(regionNames[i]) && solvers.getElectro(regionNames[i]).isSource())
        continue;
    solvers.setJToRegion(JreGlobal,regionNames[i]);
    if (solvers.isElectroHarmonic())
    {
        solvers.setJToRegion(JimGlobal,regionNames[i],true);
    }
}
solvers.setB(BreGlobal);
if (solvers.isElectroHarmonic())
{
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
