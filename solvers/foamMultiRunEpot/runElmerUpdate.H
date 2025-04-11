
// Clean global current fields before applying sources
/*if (solvers.hasElectricSources())
{
    JreGlobal = Zero;
    if (solvers.isElectroHarmonic())
    {
        JimGlobal = Zero;
    }
}*/
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
    /*if (solvers.isSource(regionNames[i]))
    {
        // For these regions current density is calculated on OpenFOAM side
        // and incorporated in Elmer as an external current source.
        const volVectorField Jre = solvers.getSourceJ(regionNames[i],false);
        solvers.vectorFieldToGlobal(
            JreGlobal,
            Jre,
            regionNames[i]
        );
        if (solvers.isElectroHarmonic())
        {
            const volVectorField Jim = solvers.getSourceJ(regionNames[i],true);
            solvers.vectorFieldToGlobal(
                JimGlobal,
                Jim,
                regionNames[i]
            );
        }
    }*/
}
//JreGlobal.write();
//JimGlobal.write();
// Send fields to Elmer
Info<< "\nSending fields to Elmer for all regions " << endl;
if (initialize_elmer) sending.initialize();
sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
sending.sendVector(UGlobal);
/*if (solvers.hasElectricSources())
{
    sending.sendVector(JreGlobal);
    if (solvers.isElectroHarmonic())
    {
        sending.sendVector(JimGlobal);
    }
}*/
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
    // Skip update if current was calculated in OpenFOAM
    solvers.setJToRegion(JreGlobal,regionNames[i]);
    if (solvers.isElectroHarmonic())
    {
        solvers.setJToRegion(JimGlobal,regionNames[i],true);
    }
    if (solvers.needsReference(regionNames[i]))
    {
        solvers.setJRefToRegion(JreGlobal,regionNames[i]);
        if (solvers.isElectroHarmonic())
        {
            solvers.setJRefToRegion(JimGlobal,regionNames[i],true);
        }
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
if (Pstream::master() && logElmerTime)
{
    std::ofstream elmerTimes(elmerTimesFileName, std::ios::app);
    if (elmerTimes.is_open())
    {
        elmerTimes << runTime.timeName() << std::endl;
        elmerTimes.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << elmerTimesFileName << " for writing!\n" << abort(FatalError);
}
// Initialization is needed only on the first update.
if (initialize_elmer)
    initialize_elmer = false;

if (!skipControllerUpdate)
    solvers.updateFeedbackControl();//JreGlobal,JimGlobal
