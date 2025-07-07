
volScalarField magU(mag(U));
scalarField volume = mesh.V();
scalar totalVolume = gSum(volume);

volScalarField Ux(mag(U & Foam::vector(1,0,0)));
volScalarField Uy(mag(U & Foam::vector(0,1,0)));
volScalarField Uz(mag(U & Foam::vector(0,0,1)));
scalar UmagIntegral = gSum(magU*volume)/totalVolume;
scalar UxIntegral = gSum(Ux*volume)/totalVolume;
scalar UyIntegral = gSum(Uy*volume)/totalVolume;
scalar UzIntegral = gSum(Uz*volume)/totalVolume;
// Write
if (Pstream::master())
{
    string velocityFileName = "postProcessing/velocityIntegrals.dat";
    std::ofstream velocityFile(velocityFileName, std::ios::app);

    if (velocityFile.is_open())
    {
        velocityFile << runTime.timeName() << "\t" << max(magU).value() << "\t" << UmagIntegral << "\t"
        << UxIntegral << "\t" << UyIntegral << "\t" << UzIntegral << std::endl;
        velocityFile.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << velocityFileName << " for writing!\n" << abort(FatalError);
}
if (regionSolver.isFluid())//fluid == hydro + thermo
{
    const volScalarField& JJsigma = regionSolver.getElectro().JJsigma;
    volScalarField T = regionSolver.getTemperature() + absoluteZero;
    scalar TIntegral = gSum(T*volume)/totalVolume;
    scalar JJsigmaIntegral = gSum(JJsigma*volume);
    if (Pstream::master())
    {
        // Write
        string thermalFileName = "postProcessing/thermalIntegrals.dat";
        std::ofstream thermalFile(thermalFileName, std::ios::app);

        if (thermalFile.is_open())
        {
            thermalFile << runTime.timeName() << "\t" << max(T).value() << "\t" << min(T).value() << "\t" << TIntegral <<
            "\t" << max(JJsigma).value() << "\t" << min(JJsigma).value() << "\t" << JJsigmaIntegral << std::endl;
            thermalFile.close();
        }
        else FatalErrorInFunction << "ERROR: Couldn't open " << thermalFileName << " for writing!\n" << abort(FatalError);
    }
}
