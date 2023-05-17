---
layout: default
title: transformToSpatialDomainWithFAllDerivatives
parent: WVTransform
grand_parent: Classes
nav_order: 172
mathjax: true
---

#  transformToSpatialDomainWithFAllDerivatives

transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.


---

## Declaration
```matlab
 [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(u_bar)
```
## Parameters
+ `u_bar`  variable with dimensions $$(k,l,j)$$

## Returns
+ `u`  variable u with dimensions $$(x,y,z)$$
+ `ux`  variable du/dx with dimensions $$(x,y,z)$$
+ `uy`  variable du/dy with dimensions $$(x,y,z)$$
+ `uz`  variable du/dz with dimensions $$(x,y,z)$$

## Discussion

This performs the same operation as `transformToSpatialDomainWithF`, but also returns the first-derivative in all three spatial directions.

The computation of these derivatives can be performed more efficiently if done simultaneously. So when performance is a requirement, it can be useful to call this function rather than request the derivatives individually.

