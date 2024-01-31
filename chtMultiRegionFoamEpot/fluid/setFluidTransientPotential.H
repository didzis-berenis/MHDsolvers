        
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
    forAll(fluidRegions, i)
    {
        fvMesh& mesh = fluidRegions[i];
        volVectorField& U = UFluid[i];
        volVectorField& U_old = U_oldFluid[i];
        volVectorField& J = JFluid[i];
        volVectorField& B = BFluid[i];
        dimensionedScalar& sigma = sigmaFluid[i];
        volScalarField& PotE = PotEFluid[i];

        volVectorField JUB = J;
        {
            #include "PotEEqn.H"
        }
        JxBFluid[i] =  ((J+JUB) ^ B );
        JJsigmaFluid[i] =  ((J+JUB) & (J+JUB))/sigma;
    }
}
