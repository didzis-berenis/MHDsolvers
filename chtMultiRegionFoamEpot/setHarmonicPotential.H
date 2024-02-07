
dimensionedScalar smallU
(
    "smallU",
    dimensionSet(0, 1, -1, 0, 0, 0 ,0),
    1e-6
);

scalar maxRemDiff_local = SMALL;        

scalar maxRelDiff_local = SMALL;

forAll(fluidRegions, i)
{
    maxRemDiff_local = max(
        magneticReynolds[i]*max(mag(U_oldFluid[i]-UFluid[i])).value(),
        maxRemDiff_local);        

    maxRelDiff_local = max(
        (max(mag(U_oldFluid[i]-UFluid[i])/(average(mag(UFluid[i]))+smallU))).value(),
        maxRemDiff_local);
}

if((maxRelDiff_local>maxRelDiff || maxRelDiff<SMALL) && maxRelDiff+SMALL<=1.0) {
    doElmer = true;
}
else if(maxRemDiff_local>maxRemDiff && maxRelDiff-SMALL<=1.0) {
    doElmer = true;
}

// Calculate electric potential if current density will not be updated
if (!doElmer)
{
    fvMesh& mesh = meshGlobal;
    volVectorField& U = UGlobal;
    volVectorField& U_old = U_oldGlobal;
    volVectorField& Jre = JreGlobal;
    volVectorField& Jim = JimGlobal;
    volVectorField& Bre = BreGlobal;
    volVectorField& Bim = BimGlobal;
    volScalarField& PotEim = PotEimGlobal;
    volScalarField& PotEre = PotEreGlobal;
    volScalarField& sigma = sigmaGlobal;

    volVectorField JUBre = Jre;
    {
        #include "PotEreEqn.H"
    }
    volVectorField JUBim = Jim;
    {
        #include "PotEimEqn.H"
    }
    JxBGlobal =  0.5*(((Jre+JUBre) ^ Bre) + ((Jim+JUBim) ^ Bim) );
	JJsigmaGlobal = 0.5*(((Jre+JUBre) & (Jre+JUBre)) + ((Jim+JUBim) & (Jim+JUBim)) )*sigmaGlobalInv;
    scalarGlobalToMultiRegion(JJsigmaGlobal,JJsigmaFluid,JJsigmaSolid);
    vectorGlobalToMultiRegion(JxBGlobal,JxBFluid,JxBSolid);
}
