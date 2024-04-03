
scalarGlobalToMultiRegion(JJsigmaGlobal,JJsigmaRegion);
vectorGlobalToMultiRegion(JxBGlobal,JxBRegion);
forAll(regionNames, i)
{
    solvers.setJxB(regionNames[i],JxBRegion[i]);
    solvers.setJJsigma(regionNames[i],JJsigmaRegion[i]);
}
