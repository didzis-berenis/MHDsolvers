
scalarField volume = meshes.completeMesh().V();
pointField coordinates = meshes.completeMesh().C();
scalar totalVolume = gSum(volume);

volVectorField curlU(fvc::curl(U));
volScalarField curlUx((curlU & Foam::vector(1,0,0)));
volScalarField curlUy((curlU & Foam::vector(0,1,0)));
volScalarField curlUz((curlU & Foam::vector(0,0,1)));

scalar curlUxIntegralAbs = gSum(mag(curlUx)*volume)/totalVolume;
scalar curlUyIntegralAbs = gSum(mag(curlUy)*volume)/totalVolume;
scalar curlUzIntegralAbs = gSum(mag(curlUz)*volume)/totalVolume;

scalar maxX, maxY, maxZ;
maxX = maxY = maxZ = -vGreat;
scalar minX, minY, minZ;
minX = minY = minZ = vGreat;
forAll (coordinates, cellI)
{
    point location = coordinates[cellI];
    scalar cellX = location.x();
    if (cellX > maxX)
    {
        maxX = cellX;
    }
    if (cellX < minX)
    {
        minX = cellX;
    }
    scalar cellY = location.y();
    if (cellY > maxY)
    {
        maxY = cellY;
    }
    if (cellY < minY)
    {
        minY = cellY;
    }
    scalar cellZ = location.z();
    if (cellZ > maxZ)
    {
        maxZ = cellZ;
    }
    if (cellZ < minZ)
    {
        minZ = cellZ;
    }
}
scalar sizeX = maxX - minX;
scalar sizeY = maxY - minY;
scalar sizeZ = maxZ - minZ;

if (curlFile.is_open())
{
    curlFile << runTimes.completeTime().timeName()
    << "\t" << max(curlUx).value()*sizeY*sizeZ << "\t" << max(curlUy).value()*sizeX*sizeZ << "\t" << max(curlUz).value()*sizeY*sizeX
    << "\t" << min(curlUx).value()*sizeY*sizeZ << "\t" << min(curlUy).value()*sizeX*sizeZ << "\t" << min(curlUz).value()*sizeY*sizeX
    << "\t" << curlUxIntegralAbs*sizeY*sizeZ << "\t" << curlUyIntegralAbs*sizeX*sizeZ << "\t" << curlUzIntegralAbs*sizeY*sizeX
    << std::endl;
}
else FatalErrorInFunction << "ERROR: Couldn't open " << curlFileName << " for writing!\n" << abort(FatalError);
