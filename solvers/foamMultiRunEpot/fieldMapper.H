//#include <map>
//#include <set>
//std::map<word, bool> keepField;

word getRegionPath(fvMesh& regionMesh)
{
    // get volScalarFields
    List<word> volScalarFieldNames(
        regionMesh.thisDb().names<volScalarField>()
        );
    HashTable<const volScalarField*> volScalarFieldTable(
        regionMesh.thisDb().lookupClass<volScalarField>()
        );
    forAll(volScalarFieldTable,i)
    {
        int writeOption = volScalarFieldTable[volScalarFieldNames[i]]->writeOpt();
        word thisRegionPath = volScalarFieldTable[volScalarFieldNames[i]]->path();
        if (writeOption == 0 && thisRegionPath != "")
        {
            //return the first valid instance
            return thisRegionPath;
        }
    }

    return word::null;
}
