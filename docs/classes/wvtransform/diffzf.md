---
layout: default
title: diffZF
parent: WVTransform
grand_parent: Classes
nav_order: 77
mathjax: true
---

#  diffZF

differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)


---

## Declaration
```matlab
 uz = diffZF(u)
```
## Parameters
+ `u`  variable with dimensions $$(x,y,z)$$

## Returns
+ `uz`  differentiated variable with dimensions $$(x,y,z)$$

## Discussion

Each subclass implements this operation differently, depending on the vertical modes being used.

For hydrostatic vertical modes with a rigid-lid and zero buoyancy anomaly,

$$
\partial_z u = -\frac{1}{g} N^2(z)  \mathcal{G}^{-1} \left[ \mathcal{F} \left[ u \right] \right]
$$

where we've used the same notation as defined for the [discrete transformations](/mathematical-introduction/transformations.html).

