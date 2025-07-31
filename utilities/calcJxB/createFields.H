// Initialize fields
//volVectorField JxB = regionSolver.getElectro().JxB;
//volScalarField JJsigma = regionSolver.getElectro().JJsigma;
volVectorField Jre = regionSolver.getElectro().J();
volVectorField Jim = regionSolver.getElectro().J(true);
volVectorField Bre = regionSolver.getElectro().B();
volVectorField Bim = regionSolver.getElectro().B(true);
const volVectorField& U = regionSolver.getU();
// Cannot use const reference for Elmer functions,
// so preparing designated field for sending velocity to Elmer.
volVectorField Usent = U;
