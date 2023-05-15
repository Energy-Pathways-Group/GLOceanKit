---
layout: default
title: transformToSpatialDomainWithGAllDerivatives
parent: WVTransform
grand_parent: Classes
nav_order: 174
mathjax: true
---

#  transformToSpatialDomainWithGAllDerivatives

transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.


---

## Declaration
```matlab
 [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives( w_bar )
```
## Parameters
+ `w_bar`  variable with dimensions $$(k,l,j)$$

## Returns
+ `w`  variable w with dimensions $$(x,y,z)$$
+ `wx`  variable dw/dx with dimensions $$(x,y,z)$$
+ `wy`  variable dw/dy with dimensions $$(x,y,z)$$
+ `wz`  variable dw/dz with dimensions $$(x,y,z)$$

## Discussion

This performs the same operation as `transformToSpatialDomainWithG`, but also returns the first-derivative in all three spatial directions.

The computation of these derivatives can be performed more efficiently if done simultaneously. So when performance is a requirement, it can be useful to call this function rather than request the derivatives individually.

