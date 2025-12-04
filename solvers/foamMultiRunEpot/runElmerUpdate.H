
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
}
for (int iter = 0; iter < nElmerCorrectors; iter++)
{
    Info << "\nElmer iteration: " << iter << endl;
    // Send fields to Elmer
    if (solvers.checkIfAnyVelocityRegions())
    {
        Info << "Sending fields to Elmer for all regions" << endl;
        if (initialize_elmer) sending.initialize();
        sending.sendStatus(elmer_status); // 1=ok, 0=lastIter, -1=error
        sending.sendVector(UGlobal);
    }
    // Receive fields from Elmer
    Info << "Receiving fields from Elmer for all regions" << endl;
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
    // Initialization is needed only on the first update.
    if (initialize_elmer)
        initialize_elmer = false;

    if (countControlSteps) controlStepsTaken++;
    //Quit the loop if bad status
    if (elmer_status != 1) break;
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
    /*if (solvers.referenceNotSet(regionNames[i]) && controlStepsTaken > 10)
    {
        solvers.setJRefToRegion(JreGlobal,regionNames[i]);
        if (solvers.isElectroHarmonic())
        {
            solvers.setJRefToRegion(JimGlobal,regionNames[i],true);
        }
    }*/
}
solvers.setB(BreGlobal);
if (solvers.isElectroHarmonic())
{
    solvers.setB(BimGlobal,true);
}
//Update region Lorentz force term and Joule heating
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
if (!skipControllerUpdate)
    solvers.updateFeedbackControl();
