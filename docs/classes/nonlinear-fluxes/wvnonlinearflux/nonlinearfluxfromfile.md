---
layout: default
title: nonlinearFluxFromFile
parent: WVNonlinearFlux
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  nonlinearFluxFromFile

initialize a nonlinear flux operation from NetCDF file


---

## Declaration
```matlab
 nlFluxOp = nonlinearFluxFromFile(ncfile,wvt)
```
## Parameters
+ `wvt`  the WVTransform to be used

## Returns
+ `nlFluxOp`  a new instance of WVNonlinearFluxOperation

## Discussion

  Subclasses to should override this method to enable model
  restarts. This method works in conjunction with -writeToFile
  to provide restart capability.
 
        
Help for WVNonlinearFlux.nonlinearFluxFromFile is inherited from superclass WVNonlinearFluxOperation