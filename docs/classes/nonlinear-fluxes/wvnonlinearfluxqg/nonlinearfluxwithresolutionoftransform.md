---
layout: default
title: nonlinearFluxWithResolutionOfTransform
parent: WVNonlinearFluxQG
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  nonlinearFluxWithResolutionOfTransform

create a new nonlinear flux operation with double the resolution


---

## Declaration
```matlab
 nlFluxOp = nonlinearFluxWithResolutionOfTransform(wvtX2)
```
## Parameters
+ `wvtX2`  the WVTransform with increased resolution

## Returns
+ `nlFluxOp`  a new instance of WVNonlinearFluxOperation

## Discussion

  Subclasses to should override this method an implement the
  correct logic. If the nonlinear flux operation makes no
  specific assumptions about the resolution of WVTransform,
  then this will work without modification.
 
        
Help for WVNonlinearFluxQG/nonlinearFluxWithResolutionOfTransform is inherited from superclass WVNonlinearFluxOperation
