#include <map>
#include <set>
std::map<word, bool> keepField;

List<Pair<word>> getFieldPaths(fvMesh& regionMesh)
{
    List<Pair<word>> fieldPaths;
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
        word fieldFileName = volScalarFieldNames[i];
        word fieldFilePath = volScalarFieldTable[volScalarFieldNames[i]]->filePath();
        if (writeOption == 0)
        {
            fieldPaths.append(Pair<word>(fieldFileName,fieldFilePath));
        }
    }
    // get volVectorFields
    List<word> volVectorFieldNames(
        regionMesh.thisDb().names<volVectorField>()
        );
    HashTable<const volVectorField*> volVectorFieldTable(
        regionMesh.thisDb().lookupClass<volVectorField>()
        );
    forAll(volVectorFieldTable,i)
    {
        int writeOption = volVectorFieldTable[volVectorFieldNames[i]]->writeOpt();
        word fieldFileName = volVectorFieldNames[i];
        word fieldFilePath = volVectorFieldTable[volVectorFieldNames[i]]->filePath();
        if (writeOption == 0)
        {
            fieldPaths.append(Pair<word>(fieldFileName,fieldFilePath));
        }
    }
    // get surfaceScalarFields
    List<word> surfaceScalarFieldNames(
        regionMesh.thisDb().names<surfaceScalarField>()
        );
    HashTable<const surfaceScalarField*> surfaceScalarFieldTable(
        regionMesh.thisDb().lookupClass<surfaceScalarField>()
        );
    forAll(surfaceScalarFieldTable,i)
    {
        int writeOption = surfaceScalarFieldTable[surfaceScalarFieldNames[i]]->writeOpt();
        word fieldFileName = surfaceScalarFieldNames[i];
        word fieldFilePath = surfaceScalarFieldTable[surfaceScalarFieldNames[i]]->filePath();
        if (writeOption == 0)
        {
            fieldPaths.append(Pair<word>(fieldFileName,fieldFilePath));
        }
    }

    return fieldPaths;
}
