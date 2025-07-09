
scalarField volume = meshes.completeMesh().V();
pointField coordinates = meshes.completeMesh().C();
scalar totalVolume = gSum(volume);

volVectorField curlU(fvc::curl(U));
volScalarField magCurlU(mag(curlU));
volScalarField curlUx((curlU & Foam::vector(1,0,0)));
volScalarField curlUy((curlU & Foam::vector(0,1,0)));
volScalarField curlUz((curlU & Foam::vector(0,0,1)));
scalar curlUxIntegral = gSum(curlUx*volume)/totalVolume;
scalar curlUyIntegral = gSum(curlUy*volume)/totalVolume;
scalar curlUzIntegral = gSum(curlUz*volume)/totalVolume;

scalar maxCurlU = -vGreat;
scalar minCurlU = vGreat;
vector maxPosition = vector(0,0,0);
vector minPosition = vector(0,0,0);
forAll (curlUx, cellI)
{
    scalar cellCurlUx = curlUx[cellI];
    scalar cellCurlUy = curlUy[cellI];
    scalar cellCurlUz = curlUz[cellI];
    scalar maxCellCurlU = max(max(cellCurlUx,cellCurlUy),cellCurlUz);
    scalar minCellCurlU = min(min(cellCurlUx,cellCurlUy),cellCurlUz);
    point location = coordinates[cellI];
    if (maxCellCurlU > maxCurlU)
    {
        maxCurlU = maxCellCurlU;
        maxPosition = location;
    }
    if (minCellCurlU < minCurlU)
    {
        minCurlU = minCellCurlU;
        minPosition = location;
    }
}
scalar signedAbsCurlU = maxCurlU;
vector signedAbsPosition = maxPosition;
int attractorState = 1;
if (abs(minCurlU)>maxCurlU)
{
    signedAbsCurlU = minCurlU;
    signedAbsPosition = minPosition;
    attractorState = 0;
}

string curlFileName = "postProcessing/curlIntegrals.dat";
std::ofstream curlFile(curlFileName, std::ios::app);

if (curlFile.is_open())
{
    curlFile << runTimes.completeTime().timeName() << "\t" << max(magCurlU).value()
    << "\t" << max(curlUx).value() << "\t" << max(curlUy).value() << "\t" << max(curlUz).value()
    << "\t" << min(curlUx).value() << "\t" << min(curlUy).value() << "\t" << min(curlUz).value()
    << "\t" << maxCurlU << "\t" << maxPosition[0]  << "\t" << maxPosition[1]  << "\t" << maxPosition[2]
    << "\t" << minCurlU << "\t" << minPosition[0]  << "\t" << minPosition[1]  << "\t" << minPosition[2]
    << "\t" << signedAbsCurlU << "\t" << signedAbsPosition[0]  << "\t" << signedAbsPosition[1]  << "\t" << signedAbsPosition[2]
    << "\t" << curlUxIntegral << "\t" << curlUyIntegral << "\t" << curlUzIntegral
    << "\t" << attractorState
    << std::endl;
    curlFile.close();
}
else FatalErrorInFunction << "ERROR: Couldn't open " << curlFileName << " for writing!\n" << abort(FatalError);
