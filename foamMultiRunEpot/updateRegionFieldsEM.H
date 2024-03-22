
scalarGlobalToMultiRegion(JJsigmaGlobal,JJsigmaRegion);
vectorGlobalToMultiRegion(JxBGlobal,JxBRegion);
forAll(regionSolverNames, i)
{
    const word& solverName = regionSolverNames[i].second();
    if (solverName == "conductingFluid")
    {
        solvers.fluids(fluidID[solverName])->setJxB(JxBRegion[i]);
        solvers.fluids(fluidID[solverName])->setJJsigma(JJsigmaRegion[i]);
    }
    if (solverName == "conductingSolid")
    {
        solvers.solids(solidID[solverName])->setJxB(JxBRegion[i]);
        solvers.solids(solidID[solverName])->setJJsigma(JJsigmaRegion[i]);
    }
}
