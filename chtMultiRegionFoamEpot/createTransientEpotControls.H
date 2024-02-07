
//Electric potential reference point
label PotERefCell = 0;
scalar PotERefValue = 0.0;
setRefCell
( 
    PotEGlobal,
    pimple.dict(),
    PotERefCell,
    PotERefValue
);
