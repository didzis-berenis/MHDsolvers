
scalarField volume = meshes.completeMesh().V();
pointField coordinates = meshes.completeMesh().C();
scalar pressureIntegral = gSum(p*volume);
scalar maxP = -vGreat;
scalar minP = vGreat;
vector maxPosition = vector(0,0,0);
vector minPosition = vector(0,0,0);
forAll (p, cellI)
{
    scalar cellPressure = p[cellI];
    point location = coordinates[cellI];
    if (cellPressure > maxP)
    {
        maxP = cellPressure;
        maxPosition = location;
    }
    if (cellPressure < minP)
    {
        minP = cellPressure;
        minPosition = location;
    }
}

//Info << "Pressure check: max " << maxP << " == " << max(p).value() << "; min "
//<< minP << " == " << min(p).value() << endl;

if (pressureFile.is_open())
{
    pressureFile << runTimes.completeTime().timeName()  << "\t" << pressureIntegral
    << "\t" << maxP << "\t" << maxPosition[0]  << "\t" << maxPosition[1]  << "\t" << maxPosition[2]
    << "\t" << minP << "\t" << minPosition[0]  << "\t" << minPosition[1]  << "\t" << minPosition[2]
    << std::endl;
}
else FatalErrorInFunction << "ERROR: Couldn't open " << pressureFileName << " for writing!\n" << abort(FatalError);
