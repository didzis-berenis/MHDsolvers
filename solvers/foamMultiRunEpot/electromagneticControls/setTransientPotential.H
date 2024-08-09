        
if (adjustableRunTime)
{
    if (runTime.writeTime()) doElmer = true;
}
else
{
    writeCounter++;
    if ( (writeCounter % writeMultiplier) == 0 && runTime.run()) doElmer = true;
}

// Calculate electric potential if current density will not be updated
if (!doElmer)
{
    solvers.setCorrectElectromagnetics();
/*
    fvMesh& mesh = meshGlobal;
    volVectorField& U = UGlobal;
    volVectorField& U_old = U_oldGlobal;
    volVectorField& J = JGlobal;
    volVectorField& B = BGlobal;
    volScalarField& PotE = PotEGlobal;
    volScalarField& sigma = sigmaGlobal;

    volVectorField JUB = J;
    {
        #include "PotEEqn.H"
    }
    JxBGlobal =  ((J+JUB) ^ B );
	JJsigmaGlobal = ((J+JUB) & (J+JUB))*sigmaGlobalInv;//multiply by inverse of sigmaGlobal
    #include "updateRegionFieldsEM.H"
    
*/
}