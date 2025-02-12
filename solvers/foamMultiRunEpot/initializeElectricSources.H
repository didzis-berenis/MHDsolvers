
Info<< nl << "Initializing electromagnetic sources\n" << endl;
bool solved = false;
solvers.setGlobalPrefix();
// Analoge for pimple.loop(), but only for selected regions.
for (int it=1; it<=pimple.nCorr(); it++)
{
    pimple.read();
    if (it==1)
        pimple.resetCorrSolveIndex();
    else
        pimple.updateCorrSolveIndex();

    forAll(regionNames, i)
    {
        if (solvers.isSource(regionNames[i]))
            pimple.pimple(solvers.getRegionId(regionNames[i])).storePrevIterFields();
    }
    forAll(regionNames, i)
    {
        if (solvers.isSource(regionNames[i]))
        {
            label thisLabel = solvers.getRegionId(regionNames[i]);
            pimple.pimple(thisLabel).storePrevIterFields();
            bool isFinal = solved || it==pimple.nCorr();
            pimple.pimple(thisLabel).updateFinal(isFinal);
        }
    }

    Info<< "Iteration " << it << endl;
    while (solvers.correctElectroPotential())
    {
        forAll(regionNames, i)
        {
            if (solvers.isSource(regionNames[i]))
                solvers.solveElectromagnetics(regionNames[i]);
        }
    }
    forAll(regionNames, i)
    {
        if (solvers.isSource(regionNames[i]))
            solvers(regionNames[i]).postCorrector();
    }
    if (solved)
    {
        solvers.setGlobalPrefix();
        Info<< "Electromagnetic source calculation converged in " << it
        << " iterations" << endl;
        break;
    }
    bool result = true;
    forAll(regionNames, i)
    {
        if (solvers.isSource(regionNames[i]))
            result = pimple.pimple(solvers.getRegionId(regionNames[i])).corrCriteriaSatisfied() && result;
    }
    solvers.setGlobalPrefix();
    solved = result;
}
if (solvers.isElectroHarmonic())
{
    // Do one extra Elmer iteration to stabilize current.
    #include "runElmerUpdate.H"
}
