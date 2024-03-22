//write selected global fields
if(runTime.writeTime())
{
    //get velocity fields from fluid domains
    forAll(regionSolverNames, i)
    {
        const word& solverName = regionSolverNames[i].second();
        if (solverName == "conductingFluid")
        {
            volVectorField& U = solvers.fluids(fluidID[solverName])->getVelocity();
            vectorFieldToGlobal(UGlobal,U,solverName);
            //get pressure fields from fluid domains
            volScalarField& p = solvers.fluids(fluidID[solverName])->getPressure();
            scalarFieldToGlobal(pGlobal,p,solverName);
            //get temperature fields from all domains
            volScalarField& T = solvers.fluids(fluidID[solverName])->getTemperature();
            scalarFieldToGlobal(TGlobal,T,solverName);
        }
        if (solverName == "conductingSolid")
        {
            volScalarField& T = solvers.solids(solidID[solverName])->getTemperature();
            scalarFieldToGlobal(TGlobal,T,solverName);
        }
    }
	UGlobal.correctBoundaryConditions();
    UGlobal.write();
	TGlobal.correctBoundaryConditions();
    TGlobal.write();
	pGlobal.correctBoundaryConditions();
    pGlobal.write();
    //joule heating and IDs should be previously assigned to global
	JJsigmaGlobal.correctBoundaryConditions();
    JJsigmaGlobal.write();
    GeometryIdsGlobal.write();

#if (ELMER_TIME == HARMONIC_TIME)
	PotEreGlobal.correctBoundaryConditions();
    PotEreGlobal.write();
	PotEimGlobal.correctBoundaryConditions();
    PotEimGlobal.write();
	JreGlobal.correctBoundaryConditions();
    JreGlobal.write();
	JimGlobal.correctBoundaryConditions();
    JimGlobal.write();
#elif (ELMER_TIME == TRANSIENT_TIME)
	PotEGlobal.correctBoundaryConditions();
    PotEGlobal.write();
	JGlobal.correctBoundaryConditions();
    JGlobal.write();
#endif
}
