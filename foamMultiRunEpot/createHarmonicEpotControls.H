
//Electric potential reference point
label PotEreRefCell = 0;
scalar PotEreRefValue = 0.0;
setRefCell
( 
    PotEreGlobal,
    pimple.dict(),
    PotEreRefCell,
    PotEreRefValue
);

label PotEimRefCell = 0;
scalar PotEimRefValue = 0.0;
setRefCell
( 
    PotEimGlobal,
    pimple.dict(),
    PotEimRefCell,
    PotEimRefValue
);
