//write selected global fields
if(runTime.writeTime())
{
    //get velocity fields from fluid domains
    forAll(regionNames, i)
    {
        const word& regionName = regionNames[i];
        if (solvers.isFluid(regionNames[i]))
        {
            URegion[i] = solvers.getVelocity(regionNames[i]);
        }
        vectorFieldToGlobal(UGlobal,URegion[i],regionName);
        //get pressure fields from fluid domains
        if (solvers.isFluid(regionNames[i]))
        {
            pRegion[i] = solvers.getPressure(regionNames[i]);
        }
        scalarFieldToGlobal(pGlobal,pRegion[i],regionName);
        //get temperature fields from all domains
        if (solvers.isFluid(regionNames[i]) || solvers.isSolid(regionNames[i]))
        {
            TRegion[i] = solvers.getTemperature(regionNames[i]);
        }
        scalarFieldToGlobal(TGlobal,TRegion[i],regionName);
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
