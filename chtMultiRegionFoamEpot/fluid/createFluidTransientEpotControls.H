
//Electric potential reference point
forAll(fluidRegions, i)
{
    label PotERefCell = 0;
    scalar PotERefValue = 0.0;
    setRefCell
    ( 
        PotEFluid[i],
        pimples.pimple(i).dict(),
        PotERefCell,
        PotERefValue
    );
}
