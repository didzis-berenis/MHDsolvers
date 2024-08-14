int fluidRegionsSize = 0;
forAll(regionNames, i)
{
    if (solvers.isFluid(regionNames[i])) fluidRegionsSize++;
}
forAll(regionNames, i)
{
    const word& regionName = regionNames[i];
    if (solvers.isFluid(regionName))
    {
        const volVectorField& U = solvers.getFluid(regionNames[i]).U;
        volScalarField magU(mag(U));
        scalarField V = solvers.mesh(regionNames[i]).V();
        scalar UmagIntegral = 0;
        scalar UxIntegral = 0;
        scalar UyIntegral = 0;
        scalar UzIntegral = 0;

        const volScalarField& T = solvers.getFluid(regionNames[i]).thermo.T();
        const volScalarField& JJsigma = solvers.getElectro(regionNames[i]).JJsigma;
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
        string regionNameStr = "";
        if (fluidRegionsSize > 1)
        {
            regionNameStr = "_" + regionName + "_";
        }
        // Write
        string velocityFileName = "postProcessing/velocityIntegrals" + regionNameStr + std::to_string(Pstream::myProcNo()) + ".dat";
        std::ofstream velocityFile(velocityFileName, std::ios::app);

        if (velocityFile.is_open())
        {
            velocityFile << runTime.timeName() << "\t" << max(magU).value() << "\t" << UmagIntegral << "\t" << UxIntegral << "\t" 
            << UyIntegral << "\t" << UzIntegral << std::endl;
            velocityFile.close();
        }
        else FatalErrorInFunction << "ERROR: Couldn't open " << velocityFileName << " for writing!\n" << abort(FatalError);

        // Write
        string thermalFileName = "postProcessing/thermalIntegrals" + regionNameStr + std::to_string(Pstream::myProcNo()) + ".dat";
        std::ofstream thermalFile(thermalFileName, std::ios::app);

        if (thermalFile.is_open())
        {
            thermalFile << runTime.timeName() << "\t" << max(T).value() << "\t" << min(T).value() << "\t" << TIntegral << 
            "\t" << max(JJsigma).value() << "\t" << min(JJsigma).value() << "\t" << JJsigmaIntegral << std::endl;
            thermalFile.close();
        }
        else FatalErrorInFunction << "ERROR: Couldn't open " << thermalFileName << " for writing!\n" << abort(FatalError);
    }

    if (solvers.isSolid(regionName))
    {
        scalarField V = solvers.mesh(regionNames[i]).V();
        const volScalarField& T = solvers.getSolid(regionNames[i]).thermo.T();
        const volScalarField& JJsigma = solvers.getElectro(regionNames[i]).JJsigma;
        scalar TIntegral = 0;
        scalar JJsigmaIntegral = 0;

        forAll (V, cellI)
        {
            scalar dV = V[cellI];
            TIntegral += T[cellI]*dV;
            JJsigmaIntegral += JJsigma[cellI]*dV;
        }
        string regionNameStr = "";
        if (regionNames.size() > 1)
        {
            regionNameStr = "_" + regionName + "_";
        }

        // Write
        string thermalFileName = "postProcessing/thermalIntegrals" + regionNameStr + std::to_string(Pstream::myProcNo()) + ".dat";
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