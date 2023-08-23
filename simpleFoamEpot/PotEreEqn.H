	//Interpolating cross product u x B over mesh faces
	surfaceScalarField psiUB = fvc::interpolate((U-U_old) ^ Bre) & mesh.Sf();
	//Poisson equation for electric potential
	fvScalarMatrix PotEEqn
	(
	fvm::laplacian(PotEre)
	==
	fvc::div(psiUB)
	);
	//Reference potential
	PotEEqn.setReference(PotEreRefCell, PotEreRefValue);
	//Solving Poisson equation
	PotEEqn.solve();

	//Computation of current density at cell faces
	surfaceScalarField En = -(fvc::snGrad(PotEre) * mesh.magSf()) + psiUB;
	//Current density at face center
	surfaceVectorField Env = En * mesh.Cf();
	//Interpolation of current density at cell center
	JUBre = sigma*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh.C()) );
	//Update current density distribution and boundary condition
	JUBre.correctBoundaryConditions();
