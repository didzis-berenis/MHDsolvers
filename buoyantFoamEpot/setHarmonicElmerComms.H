
receiving.recvVector(Jre);
receiving.recvVector(Jim);
receiving.recvVector(Bre);
receiving.recvVector(Bim);

//Lorentz force term initialization
JxB =  0.5*((Jre ^ Bre) + (Jim ^ Bim) );
JJsigma =  0.5*((Jre & Jre) + (Jim & Jim) )/sigma;
