#include <map>

std::map<std::pair<word, label>, label> localToGlobalID;
//std::map<word, label> regionID;
List<word> regionNames;
//Assigns fluid and solid region values from global to each region field
void vectorGlobalToMultiRegion(volVectorField& global,PtrList<volVectorField>& regions)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            regions[i][cellI] = global[localToGlobalID[std::make_pair(regionNames[i],cellI)]];
        }
    }
}

void scalarGlobalToMultiRegion(volScalarField& global,PtrList<volScalarField>& regions)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            regions[i][cellI] = global[localToGlobalID[std::make_pair(regionNames[i],cellI)]];
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
            global[localToGlobalID[std::make_pair(regionNames[i],cellI)]] = regions[i][cellI];
        }
    }
}

void scalarMultiRegionToGlobal(volScalarField& global,PtrList<volScalarField>& regions)
{
    forAll(regionNames, i)
    {
        forAll(regions[i], cellI)
        {
            global[localToGlobalID[std::make_pair(regionNames[i],cellI)]] = regions[i][cellI];
        }
    }
}
//Assigns single fluid/solid region values from region to global field

void scalarFieldToGlobal(volScalarField& global,volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[localToGlobalID[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}

void vectorFieldToGlobal(volVectorField& global,volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        global[localToGlobalID[std::make_pair(regionName,cellI)]] = region[cellI];
    }
}
//Assigns single fluid/solid region values from global field to region
void scalarGlobalToField(volScalarField& global,volScalarField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[localToGlobalID[std::make_pair(regionName,cellI)]];
    }
}

void vectorGlobalToField(volVectorField& global,volVectorField& region,const word& regionName)
{
    forAll(region, cellI)
    {
        region[cellI] = global[localToGlobalID[std::make_pair(regionName,cellI)]];
    }
}
