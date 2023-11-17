## v0.1: Initially included solvers. ##
* buoyantFoamMHD
* pimpleFoamMHD
* simpleFoamMHD
## v0.2: Added simpleFoam solver that updates electrical potential w/o Elmer. ##
* simpleFoamEpot
## v0.2.1: Added pimpleFoam and buoyantFoam solvers that update electrical potential w/o Elmer. ##
* buoyantFoamEpot
* pimpleFoamEpot
## v0.2.2: Added generateSimpleEOFFiles and generatePimpleEOFFiles for populating case directory with O2E files. ##
* generateSimpleEOFFiles
* generatePimpleEOFFiles
## v0.2.3: Added pimpleFoam and buoyantFoam solvers that updates through Elmer after fixed time interval, but updates only electrical potential for the time steps in-between. ##
* buoyantFoamEpotTransient
* pimpleFoamEpotTransient
## v0.2.4: Fixed generateSimpleEOFFiles. ##
## v0.2.5: Added initializeTemperature for generating initial temperature gradient ##
* initializeTemperature
