---
layout: default
title: diffY
parent: WVTransform
grand_parent: Classes
nav_order: 82
mathjax: true
---

#  diffY

differentiate a spatial variable in the y-direction


---

## Declaration
```matlab
 du = diffY(u,n)
```
## Parameters
+ `u`  variable with dimensions $$(x,y,z)$$
+ `n`  (optional) order of differentiation $$\frac{d^n}{dy^n}$$ (default 1)

## Returns
+ `du`  differentiated variable in the spatial domain

## Discussion

  Performs spectral differentiation on variable u.
 
          
