
scalarField volume = meshes.completeMesh().V();
scalar totalVolume = gSum(volume);

volVectorField curlU(fvc::curl(U));
volScalarField curlUx((curlU & Foam::vector(1,0,0)));
volScalarField curlUy((curlU & Foam::vector(0,1,0)));
volScalarField curlUz((curlU & Foam::vector(0,0,1)));
scalar curlUxIntegral = gSum(curlUx*volume)/totalVolume;
scalar curlUyIntegral = gSum(curlUy*volume)/totalVolume;
scalar curlUzIntegral = gSum(curlUz*volume)/totalVolume;

string curlFileName = "postProcessing/curlIntegrals.dat";
std::ofstream curlFile(curlFileName, std::ios::app);

if (curlFile.is_open())
{
    curlFile << runTimes.completeTime().timeName()
    << "\t" << curlUxIntegral << "\t" << curlUyIntegral << "\t" << curlUzIntegral
    << std::endl;
    curlFile.close();
}
else FatalErrorInFunction << "ERROR: Couldn't open " << curlFileName << " for writing!\n" << abort(FatalError);
