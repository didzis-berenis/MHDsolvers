PtrList<pressureReference> pressureReferenceFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    pressureReferenceFluid.set
    (
        i,
        new pressureReference
        (
            p_rghFluid[i],
            pimples.pimple(i).dict(),
            false
        )
    );

    hydrostaticInitialisation
    (
        p_rghFluid[i],
        thermoFluid[i].p(),
        rhoFluid[i],
        UFluid[i],
        ghFluid[i],
        ghfFluid[i],
        pRefFluid[i],
        thermoFluid[i],
        pimples.pimple(i).dict()
    );

    //Electric potential reference point
    label PotEreRefCell = 0;
    scalar PotEreRefValue = 0.0;
    setRefCell
    ( 
        PotEreFluid[i],
        pimples.pimple(i).dict(),
        PotEreRefCell,
        PotEreRefValue
    );
    
    label PotEimRefCell = 0;
    scalar PotEimRefValue = 0.0;
    setRefCell
    ( 
        PotEimFluid[i],
        pimples.pimple(i).dict(),
        PotEimRefCell,
        PotEimRefValue
    );
}
