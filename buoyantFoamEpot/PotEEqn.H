	//Interpolating cross product u x B over mesh faces
	surfaceScalarField psiUB = fvc::interpolate((U-U_old) ^ B) & mesh.Sf();
	//Poisson equation for electric potential
	fvScalarMatrix PotEEqn
	(
	fvm::laplacian(PotE)
	==
	fvc::div(psiUB)
	);
	//Reference potential
	PotEEqn.setReference(PotERefCell, PotERefValue);
	//Solving Poisson equation
	PotEEqn.solve();

	//Computation of current density at cell faces
	surfaceScalarField En = -(fvc::snGrad(PotE) * mesh.magSf()) + psiUB;
	//Current density at face center
	surfaceVectorField Env = En * mesh.Cf();
	//Interpolation of current density at cell center
	JUB = sigma*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh.C()) );
	//Update current density distribution and boundary condition
	JUB.correctBoundaryConditions();
