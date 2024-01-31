	//Interpolating cross product u x B over mesh faces
	surfaceScalarField psiUB = fvc::interpolate((U-U_old) ^ Bim) & mesh.Sf();
	//Poisson equation for electric potential
	fvScalarMatrix PotEEqn
	(
	fvm::laplacian(PotEim)
	==
	fvc::div(psiUB)
	);
	//Reference potential
    label PotEimRefCell = 0;
    scalar PotEimRefValue = 0.0;
	PotEEqn.setReference(PotEimRefCell, PotEimRefValue);
	//Solving Poisson equation
	PotEEqn.solve();

	//Computation of current density at cell faces
	surfaceScalarField En = -(fvc::snGrad(PotEim) * mesh.magSf()) + psiUB;
	//Current density at face center
	surfaceVectorField Env = En * mesh.Cf();
	//Interpolation of current density at cell center
	JUBim = sigma*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh.C()) );
	//Update current density distribution and boundary condition
	JUBim.correctBoundaryConditions();
