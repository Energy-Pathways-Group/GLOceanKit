---
layout: default
title: QGPVE
parent: QGPVE
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  QGPVE

initialize 3D quasigeostrophic potential vorticity flux


---

## Declaration
```matlab
 nlFlux = QGPVE(wvt,options)
```
## Parameters
+ `wvt`  a WVTransform instance
+ `shouldUseBeta`  (optional) a Boolean indicating whether or not to include beta in the flux
+ `u_damp`  (optional) characteristic speed used to set the damping. Try using wvt.uMax
+ `r`  (optional) bottom friction
+ `nu`  (optional) coefficient for damping

## Returns
+ `nlFlux`  a QGPVE instance

## Discussion

              
