---
layout: default
title: WVNonlinearFluxQGForced
parent: WVNonlinearFluxQGForced
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  WVNonlinearFluxQGForced

initialize 3D quasigeostrophic potential vorticity flux


---

## Declaration
```matlab
 nlFlux = WVNonlinearFluxQGForced(wvt,options)
```
## Parameters
+ `wvt`  a WVTransform instance
+ `shouldUseBeta`  (optional) a Boolean indicating whether or not to include beta in the flux
+ `uv_damp`  (optional) characteristic speed used to set the damping. Try using wvt.uMax
+ `r`  (optional) bottom friction
+ `nu_xy`  (optional) coefficient for damping

## Returns
+ `nlFlux`  a QGPVE instance

## Discussion

              
