const fvMesh& mesh = solidRegions[i];

solidThermo& thermo = thermoSolid[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();
volScalarField& JJsigma = JJsigmaSolid[i];

volScalarField& e = thermo.he();

const Foam::fvModels& fvModels = fvModelsSolid[i];
Foam::fvConstraints& fvConstraints = fvConstraintsSolid[i];

#include "checkRadiationModel.H"

solidNoLoopControl& pimple = pimples.solid(i);
