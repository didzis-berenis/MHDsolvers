//Assigns fluid and solid region field values from global to each region mesh
void vectorGlobalToMultiRegion(volVectorField& global,PtrList<volVectorField>& fluids,PtrList<volVectorField>& solids)
{
    int cellGlobal = 0;
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            fluids[i][cellI] = global[cellGlobal];
            cellGlobal++;
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            solids[i][cellI] = global[cellGlobal];
            cellGlobal++;
        }
    }
}

void scalarGlobalToMultiRegion(volScalarField& global,PtrList<volScalarField>& fluids,PtrList<volScalarField>& solids)
{
    int cellGlobal = 0;
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            fluids[i][cellI] = global[cellGlobal];
            cellGlobal++;
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            solids[i][cellI] = global[cellGlobal];
            cellGlobal++;
        }
    }
}
//Assigns fluid and solid region field values from each region to global mesh
void vectorMultiRegionToGlobal(volVectorField& global,PtrList<volVectorField>& fluids,PtrList<volVectorField>& solids)
{
    int cellGlobal = 0;
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            global[cellGlobal] = fluids[i][cellI];
            cellGlobal++;
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            global[cellGlobal] = solids[i][cellI];
            cellGlobal++;
        }
    }
}

void scalarMultiRegionToGlobal(volScalarField& global,PtrList<volScalarField>& fluids,PtrList<volScalarField>& solids)
{
    int cellGlobal = 0;
    forAll(fluids, i)
    {
        forAll(fluids[i], cellI)
        {
            global[cellGlobal] = fluids[i][cellI];
            cellGlobal++;
        }
    }
    forAll(solids, i)
    {
        forAll(solids[i], cellI)
        {
            global[cellGlobal] = solids[i][cellI];
            cellGlobal++;
        }
    }
}
