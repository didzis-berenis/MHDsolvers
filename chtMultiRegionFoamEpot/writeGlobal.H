//write selected global fields
if(runTime.writeTime())
{
    //get velocity fields from fluid domains
    forAll(fluidRegions, i)
    {
        vectorFluidFieldToGlobal(UGlobal,UFluid[i],i);
    }
	UGlobal.correctBoundaryConditions();
    UGlobal.write();
    //get temperature and pressure fields from fluid domains
    forAll(fluidRegions, i)
    {
        fluidThermo& thermo = thermoFluid[i];
        scalarFluidFieldToGlobal(TGlobal,thermo.T(),i);
        scalarFluidFieldToGlobal(pGlobal,thermo.p(),i);
    }
    //get temperature fields from fluid domains
    forAll(solidRegions, i)
    {
        solidThermo& thermo = thermoSolid[i];
        scalarSolidFieldToGlobal(TGlobal,thermo.T(),i);
    }
	TGlobal.correctBoundaryConditions();
    TGlobal.write();
	pGlobal.correctBoundaryConditions();
    pGlobal.write();
    //joule heating and IDs should be previously assigned to global
	JJsigmaGlobal.correctBoundaryConditions();
    JJsigmaGlobal.write();
    GeometryIdsGlobal.write();
}
