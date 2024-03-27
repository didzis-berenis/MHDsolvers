
if (adjustableRunTime)
{
    if (runTime.writeTime()) doElmer = true;
}
else
{
    writeCounter++;
    if ( (writeCounter % writeMultiplier) == 0 && runTime.run()) doElmer = true;
}

// Calculate electric potential if current density will not be updated
if (!doElmer)
{
    volVectorField JUB = J;
    {
        #include "PotEEqn.H"
    }
    JxB =  ((J+JUB) ^ B);
    regionSolver.setJxB(JxB);
    if (regionSolver.isFluid())
    {
        JJsigma =  ((J+JUB) & (J+JUB))/sigma;
        regionSolver.setJJsigma(JJsigma);
    }
}
