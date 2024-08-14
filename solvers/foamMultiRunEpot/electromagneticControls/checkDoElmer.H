
if (solvers.isElectroHarmonic())
{
    dimensionedScalar smallU
    (
        "smallU",
        dimensionSet(0, 1, -1, 0, 0, 0 ,0),
        1e-6
    );

    scalar maxRemDiff_local = SMALL;        

    scalar maxRelDiff_local = SMALL;

    forAll(regionNames, i)
    {
        if (solvers.isFluid(regionNames[i]))
        {
            URegion[i] = solvers.getFluid(regionNames[i]).U;
        }
        volVectorField& U = URegion[i];
        maxRemDiff_local = max(
            magneticReynoldsRegion[i]*max(mag(U_oldRegion[i]-U)).value(),
            maxRemDiff_local);        

        maxRelDiff_local = max(
            (max(mag(U_oldRegion[i]-U)/(average(mag(U))+smallU))).value(),
            maxRemDiff_local);
    }

    if((maxRelDiff_local>maxRelDiff || maxRelDiff<SMALL) && maxRelDiff+SMALL<=1.0) {
        doElmer = true;
    }
    else if(maxRemDiff_local>maxRemDiff && maxRelDiff-SMALL<=1.0) {
        doElmer = true;
    }
}
else
{
    if (adjustableRunTime)
    {
        if (runTime.writeTime()) doElmer = true;
    }
    else
    {
        writeCounter++;
        if ( (writeCounter % writeMultiplier) == 0 && runTime.run()) doElmer = true;
    }
}
