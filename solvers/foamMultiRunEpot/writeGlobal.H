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
    if (!solidificationGlobalEnabled)
    {
	    TGlobal.correctBoundaryConditions();
    }
    TGlobal.write();
	pGlobal.correctBoundaryConditions();
    pGlobal.write();
    //joule heating and IDs should be previously assigned to global
	JJsigmaGlobal.correctBoundaryConditions();
    JJsigmaGlobal.write();
	JxBGlobal.correctBoundaryConditions();
    JxBGlobal.write();
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

    // Cleanup
    forAll(regionNames, i)
    {
        //Pout << "Deleting folder " << regionPaths[i] << " was " <<
        fileHandler().rmDir(regionPaths[i]);//<< endl;
    }
    forAll(fieldPaths, i)
    {
        if (!keepField[fieldPaths[i].first()])
        {
            //Pout << "Deleting file " << fieldPaths[i].first() <<" was " <<
            fileHandler().rm(fieldPaths[i].second());//<< endl;
        }
    }
}
