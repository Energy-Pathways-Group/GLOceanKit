---
layout: default
title: nonlinearFluxFromFile
parent: WVNonlinearFluxQGForced
grand_parent: Classes
nav_order: 6
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
 
        
Help for WVNonlinearFluxQGForced.nonlinearFluxFromFile is inherited from superclass WVNonlinearFluxOperation
