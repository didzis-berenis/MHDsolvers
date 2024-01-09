PtrList<pressureReference> pressureReferenceFluid(fluidRegions.size());

//Electric potential reference point
forAll(fluidRegions, i)
{
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
