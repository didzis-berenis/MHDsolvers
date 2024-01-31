# About
* MHD solvers are modified OpenFOAM v10 solvers (See [OpenFOAM GitHub repository](https://github.com/OpenFOAM/OpenFOAM-10)).
* The solvers are to be used with [EOF-Library](https://github.com/jvencels/EOF-Library). Prerequisites are [OpenFOAM v10](https://openfoam.org/version/10/) and [Elmer FEM](https://www.csc.fi/web/elmer) software.
* To install all solvers, just execute ./Allwmake in a bash terminal of a properly prepared Linux environment.

## chtMultiRegionFoamEpot ##
* A modified version of OpenFOAM v10 solver chtMultiRegionFoam.
* Solver for steady or transient electromagnetically forced fluid flow and solid heat conduction, with conjugate heat transfer between regions, buoyancy effects, turbulence, reactions and radiation modelling.

## buoyantFoamEpot ##
* A modified version of OpenFOAM v10 solver buoyantFoam.
* Solver for steady or transient buoyant, turbulent flow of compressible fluids for electromagnetically forced and heated flows.

## pimpleFoamEpot ##
* A modified version of OpenFOAM v10 solver pimpleFoam.
* Transient solver for electromagnetically forced incompressible, turbulent flow of Newtonian fluids.

## simpleFoamEpot ##
* A modified version of OpenFOAM v10 solver simpleFoam.
* Steady-state solver for electromagnetically forced incompressible, turbulent flow.

