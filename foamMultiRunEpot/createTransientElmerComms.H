// Get control parameters
int writeMultiplier( readScalar(runTime.controlDict().lookup("writeMultiplier")) );
long long int writeCounter = 0;
word writeControlDict(runTime.controlDict().lookup("writeControl"));
const bool adjustableRunTime = (writeControlDict=="adjustableRunTime");
