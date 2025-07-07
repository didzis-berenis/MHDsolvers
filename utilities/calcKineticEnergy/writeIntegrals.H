
volScalarField energy((U & U));
scalarField volume = meshes.completeMesh().V();
scalar totalVolume = gSum(volume);
scalar energyIntegral = 0.5*gSum(energy*volume)/totalVolume;

string energyFileName = "postProcessing/energyIntegrals.dat";
std::ofstream energyFile(energyFileName, std::ios::app);

if (energyFile.is_open())
{
    energyFile << runTimes.completeTime().timeName() << "\t" << energyIntegral
    << std::endl;
    energyFile.close();
}
else FatalErrorInFunction << "ERROR: Couldn't open " << energyFileName << " for writing!\n" << abort(FatalError);
