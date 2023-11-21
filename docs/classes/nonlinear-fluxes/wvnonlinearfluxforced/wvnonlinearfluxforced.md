---
layout: default
title: WVNonlinearFluxForced
parent: WVNonlinearFluxForced
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  WVNonlinearFluxForced

initialize WVNonlinearFluxForced


---

## Declaration
```matlab
 nlFlux = WVNonlinearFluxForced(wvt,options)
```
## Parameters
+ `wvt`  a WVTransform instance
+ `uv_damp`  (optional) characteristic speed used to set the damping. Try using wvt.uMax.
+ `w_damp`  (optional) characteristic speed used to set the damping. Try using wvt.wMax.
+ `nu_xy`  (optional) coefficient for damping
+ `nu_z`  (optional) coefficient for damping
+ `shouldAntialias`  (optional) a Boolean indicating whether or not to antialias (default 1)

## Returns
+ `nlFlux`  a WVNonlinearFlux instance

## Discussion

                  
