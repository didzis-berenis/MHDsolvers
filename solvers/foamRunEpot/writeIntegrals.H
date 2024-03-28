
volScalarField magU(mag(U));
scalarField V = mesh.V();
scalar UmagIntegral = 0;
scalar UxIntegral = 0;
scalar UyIntegral = 0;
scalar UzIntegral = 0;

forAll (V, cellI)
{
    scalar dV = V[cellI];
    UmagIntegral += magU[cellI]*dV;
    UxIntegral += std::abs(U[cellI][0])*dV;
    UyIntegral += std::abs(U[cellI][1])*dV;
    UzIntegral += std::abs(U[cellI][2])*dV;
}
// Write
string velocityFileName = "postProcessing/velocityIntegrals" + std::to_string(Pstream::myProcNo()) + ".dat";
std::ofstream velocityFile(velocityFileName, std::ios::app);

if (velocityFile.is_open())
{
    velocityFile << runTime.timeName() << "\t" << max(magU).value() << "\t" << UmagIntegral << "\t" << UxIntegral << "\t" 
    << UyIntegral << "\t" << UzIntegral << std::endl;
    velocityFile.close();
}
else FatalErrorInFunction << "ERROR: Couldn't open " << velocityFileName << " for writing!\n" << abort(FatalError);

if (regionSolver.isFluid())
{
    volScalarField& T = regionSolver.getTemperature();
    scalar TIntegral = 0;
    scalar JJsigmaIntegral = 0;

    forAll (V, cellI)
    {
        scalar dV = V[cellI];

        TIntegral += T[cellI]*dV;
        JJsigmaIntegral += JJsigma[cellI]*dV;
    }
    // Write
    string thermalFileName = "postProcessing/thermalIntegrals" + std::to_string(Pstream::myProcNo()) + ".dat";
    std::ofstream thermalFile(thermalFileName, std::ios::app);

    if (thermalFile.is_open())
    {
        thermalFile << runTime.timeName() << "\t" << max(T).value() << "\t" << min(T).value() << "\t" << TIntegral << 
        "\t" << max(JJsigma).value() << "\t" << min(JJsigma).value() << "\t" << JJsigmaIntegral << std::endl;
        thermalFile.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << thermalFileName << " for writing!\n" << abort(FatalError);
}
