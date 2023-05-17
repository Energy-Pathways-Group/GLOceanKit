---
layout: default
title: Boussinesq
parent: Boussinesq
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  Boussinesq

initialize the Boussinesq nonlinear flux


---

## Declaration
```matlab
 nlFlux = Boussinesq(wvt,options)
```
## Parameters
+ `wvt`  a WVTransform instance
+ `uv_damp`  (optional) characteristic speed used to set the damping. Try using wvt.uMax.
+ `w_damp`  (optional) characteristic speed used to set the damping. Try using wvt.wMax.
+ `nu_xy`  (optional) coefficient for damping
+ `nu_z`  (optional) coefficient for damping
+ `shouldAntialias`  (optional) a Boolean indicating whether or not to antialias (default 1)

## Returns
+ `nlFlux`  a Boussinesq instance

## Discussion

                
