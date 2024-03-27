#include <map>

std::map<word, label> regionOffset;
//std::map<word, label> regionID;
List<word> regionNames;
//Assigns fluid and solid region values from global to each region field
void vectorGlobalToMultiRegion(volVectorField& global,PtrList<volVectorField>& regions)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            regions[i][cellI] = global[cellI+regionOffset[regionNames[i]]];
        }
    }
}

void scalarGlobalToMultiRegion(volScalarField& global,PtrList<volScalarField>& regions)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            regions[i][cellI] = global[cellI+regionOffset[regionNames[i]]];
        }
    }
}
//Assigns fluid and solid region values from each region to global field
void vectorMultiRegionToGlobal(volVectorField& global,PtrList<volVectorField>& regions)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            global[cellI+regionOffset[regionNames[i]]] = regions[i][cellI];
        }
    }
}

void scalarMultiRegionToGlobal(volScalarField& global,PtrList<volScalarField>& regions)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            global[cellI+regionOffset[regionNames[i]]] = regions[i][cellI];
        }
    }
}
//Assigns single fluid/solid region values from region to global field

void scalarFieldToGlobal(volScalarField& global,volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[cellI+regionOffset[regionName]] = region[cellI];
    }
}

void vectorFieldToGlobal(volVectorField& global,volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[cellI+regionOffset[regionName]] = region[cellI];
    }
}
//Assigns single fluid/solid region values from global field to region
void scalarGlobalToField(volScalarField& global,volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[cellI+regionOffset[regionName]];
    }
}

void vectorGlobalToField(volVectorField& global,volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[cellI+regionOffset[regionName]];
    }
}
