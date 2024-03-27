
receiving.recvVector(J);
receiving.recvVector(B);

//Lorentz force term initialization
JxB =  (J ^ B);
regionSolver.setJxB(JxB);
if (regionSolver.isFluid())
{
    JJsigma =  (J & J)/sigma;
    regionSolver.setJJsigma(JJsigma);
}
