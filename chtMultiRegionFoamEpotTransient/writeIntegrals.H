
forAll(fluidRegions, i)
{
    volVectorField& U = UFluid[i];
    volScalarField magU(mag(U));
    scalarField V = fluidRegions[i].V();
    scalar UmagIntegral = 0;
    scalar UxIntegral = 0;
    scalar UyIntegral = 0;
    scalar UzIntegral = 0;

    volScalarField& T = thermoFluid[i].T();
    volScalarField& JJsigma = JJsigmaFluid[i];
    scalar TIntegral = 0;
    scalar JJsigmaIntegral = 0;

    forAll (V, cellI)
    {
        scalar dV = V[cellI];
        UmagIntegral += magU[cellI]*dV;
        UxIntegral += std::abs(U[cellI][0])*dV;
        UyIntegral += std::abs(U[cellI][1])*dV;
        UzIntegral += std::abs(U[cellI][2])*dV;

        TIntegral += T[cellI]*dV;
        JJsigmaIntegral += JJsigma[cellI]*dV;
    }
    string regionName = "";
    if (fluidRegions.size() + i > 1)
    {
        regionName = fluidNames[i];
    }
    // Write
    string velocityFileName = "velocityIntegrals" + std::to_string(Pstream::myProcNo()) + regionName + ".dat";
    std::ofstream velocityFile(velocityFileName, std::ios::app);

    if (velocityFile.is_open())
    {
        velocityFile << runTime.timeName() << "\t" << max(magU).value() << "\t" << UmagIntegral << "\t" << UxIntegral << "\t" 
        << UyIntegral << "\t" << UzIntegral << std::endl;
        velocityFile.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << velocityFileName << " for writing!\n" << abort(FatalError);

    // Write
    string thermalFileName = "thermalIntegrals" + regionName + std::to_string(Pstream::myProcNo()) + ".dat";
    std::ofstream thermalFile(thermalFileName, std::ios::app);

    if (thermalFile.is_open())
    {
        thermalFile << runTime.timeName() << "\t" << max(T).value() << "\t" << min(T).value() << "\t" << TIntegral << 
        "\t" << max(JJsigma).value() << "\t" << min(JJsigma).value() << "\t" << JJsigmaIntegral << std::endl;
        thermalFile.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << thermalFileName << " for writing!\n" << abort(FatalError);
}

forAll(solidRegions, i)
{
    scalarField V = solidRegions[i].V();
    volScalarField& T = thermoSolid[i].T();
    volScalarField& JJsigma = JJsigmaSolid[i];
    scalar TIntegral = 0;
    scalar JJsigmaIntegral = 0;

    forAll (V, cellI)
    {
        scalar dV = V[cellI];
        TIntegral += T[cellI]*dV;
        JJsigmaIntegral += JJsigma[cellI]*dV;
    }
    string regionName = "";
    if (fluidRegions.size() + solidRegions.size() + i > 1)
    {
        regionName = solidNames[i];
    }

    // Write
    string thermalFileName = "thermalIntegrals" + regionName + std::to_string(Pstream::myProcNo()) + ".dat";
    std::ofstream thermalFile(thermalFileName, std::ios::app);

    if (thermalFile.is_open())
    {
        thermalFile << runTime.timeName() << "\t" << max(T).value() << "\t" << min(T).value() << "\t" << TIntegral << 
        "\t" << max(JJsigma).value() << "\t" << min(JJsigma).value() << "\t" << JJsigmaIntegral << std::endl;
        thermalFile.close();
    }
    else FatalErrorInFunction << "ERROR: Couldn't open " << thermalFileName << " for writing!\n" << abort(FatalError);
}
