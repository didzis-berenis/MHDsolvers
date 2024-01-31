
dimensionedScalar smallU
(
    "smallU",
    dimensionSet(0, 1, -1, 0, 0, 0 ,0),
    1e-6
);

scalar maxRemDiff_local = Rem0*max(mag(U_old-U)).value();        

scalar maxRelDiff_local = (max(mag(U_old-U)/(average(mag(U))+smallU))).value();

if((maxRelDiff_local>maxRelDiff || maxRelDiff<SMALL) && maxRelDiff+SMALL<=1.0) {
    doElmer = true;
}
else if(maxRemDiff_local>maxRemDiff && maxRelDiff-SMALL<=1.0) {
    doElmer = true;
}

// Calculate electric potential if current density will not be updated
if (!doElmer)
{
    volVectorField JUBre = Jre;
    {
        #include "PotEreEqn.H"
    }
    volVectorField JUBim = Jim;
    {
        #include "PotEimEqn.H"
    }
    JxB =  0.5*(((Jre+JUBre) ^ Bre) + ((Jim+JUBim) ^ Bim) );
}
