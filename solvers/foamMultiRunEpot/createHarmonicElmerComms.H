// Get control parameters

// Maximum allowable magnetic Reynolds number difference comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     (0,inf) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRemDiff(readScalar(runTime.controlDict().lookup("maxRemDiff")));

//
// Maximum allowable relative field difference in any cell comparing
// to last Elmer update. This option controls frequency Elmer is called.
//     >1  - once
//     1 - magnetic Reynolds number is used instead
//     (0,1) - when relative difference in any cell exceeds given value
//     0     - every iteration

scalar maxRelDiff(readScalar(runTime.controlDict().lookup("maxRelDiff")));
