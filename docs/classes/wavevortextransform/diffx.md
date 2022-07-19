---
layout: default
title: diffX
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 80
mathjax: true
---

#  diffX

differentiate a spatial variable in the x-direction


---

## Declaration
```matlab
 du = diffX(u,n)
```
## Parameters
+ `u`  variable with dimensions $$(x,y,z)$$
+ `n`  (optional) order of differentiation d^n/dx^n (default 1)

## Returns
+ `du`  differentiated variable in the spatial domain

## Discussion

  Performs spectral differentiation on variable u.
 
          
