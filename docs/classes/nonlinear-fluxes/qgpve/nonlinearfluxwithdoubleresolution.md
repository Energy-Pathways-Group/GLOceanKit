---
layout: default
title: nonlinearFluxWithDoubleResolution
parent: QGPVE
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  nonlinearFluxWithDoubleResolution

create a new nonlinear flux operation with double the resolution


---

## Declaration
```matlab
 nlFluxOp = nonlinearFluxWithDoubleResolution(wvtX2)
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
 
        
Help for QGPVE/nonlinearFluxWithDoubleResolution is inherited from superclass WVNonlinearFluxOperation