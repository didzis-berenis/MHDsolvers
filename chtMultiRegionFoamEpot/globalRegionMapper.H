#include <map>
// maps fluid region ID to cell offset in the global mesh
std::map<label, label> fluidOffset;
// maps solid region ID to cell offset in the global mesh
std::map<label, label> solidOffset;
//Assigns fluid and solid region values from global to each region field
void vectorGlobalToMultiRegion(volVectorField& global,PtrList<volVectorField>& fluids,PtrList<volVectorField>& solids)
{
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            fluids[i][cellI] = global[cellI+fluidOffset[i]];
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            solids[i][cellI] = global[cellI+solidOffset[i]];
        }
    }
}

void scalarGlobalToMultiRegion(volScalarField& global,PtrList<volScalarField>& fluids,PtrList<volScalarField>& solids)
{
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            fluids[i][cellI] = global[cellI+fluidOffset[i]];
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            solids[i][cellI] = global[cellI+solidOffset[i]];
        }
    }
}
//Assigns fluid and solid region values from each region to global field
void vectorMultiRegionToGlobal(volVectorField& global,PtrList<volVectorField>& fluids,PtrList<volVectorField>& solids)
{
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            global[cellI+fluidOffset[i]] = fluids[i][cellI];
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            global[cellI+solidOffset[i]] = solids[i][cellI];
        }
    }
}

void scalarMultiRegionToGlobal(volScalarField& global,PtrList<volScalarField>& fluids,PtrList<volScalarField>& solids)
{
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            global[cellI+fluidOffset[i]] = fluids[i][cellI];
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            global[cellI+solidOffset[i]] = solids[i][cellI];
        }
    }
}
//Assigns single fluid/solid region values from region to global field
void scalarFluidFieldToGlobal(volScalarField& global,volScalarField& fluid,label i)
{
    forAll(fluid, cellI)
    {
        global[cellI+fluidOffset[i]] = fluid[cellI];
    }
}

void scalarSolidFieldToGlobal(volScalarField& global,volScalarField& solid,label i)
{
    forAll(solid, cellI)
    {
        global[cellI+solidOffset[i]] = solid[cellI];
    }
}

void vectorFluidFieldToGlobal(volVectorField& global,volVectorField& fluid,label i)
{
    forAll(fluid, cellI)
    {
        global[cellI+fluidOffset[i]] = fluid[cellI];
    }
}

void vectorSolidFieldToGlobal(volVectorField& global,volVectorField& solid,label i)
{
    forAll(solid, cellI)
    {
        global[cellI+solidOffset[i]] = solid[cellI];
    }
}
