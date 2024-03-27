
regionSolver.getVelocity(U);
regionSolver.getPressure(p);
if (regionSolver.isFluid())
{
    regionSolver.getTemperature(T);
}
