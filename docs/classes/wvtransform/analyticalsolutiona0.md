---
layout: default
title: analyticalSolutionA0
parent: WVTransform
grand_parent: Classes
nav_order: 77
mathjax: true
---

#  analyticalSolutionA0

return the analytical solution of the geostrophic mode


---

## Declaration
```matlab
 [k,l] = analyticalSolutionA0(self)
```
## Parameters
+ `kMode`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
+ `A0`  (optional) amplitude in m. Default will use the current value.
+ `phi`  (optional) phase in radians, (0 <= phi <= 2*pi)
+ `eta`  (optional) isopycnal displacement (m)

## Returns
+ `u`  fluid velocity, u = @(x,y,z,t)
+ `v`  fluid velocity, v = @(x,y,z,t)
+ `w`  fluid velocity, w = @(x,y,z,t)
+ `eta`  isopycnal displacement, eta = @(x,y,z,t)
+ `p`  pressure, p = @(x,y,z,t)

## Discussion

  Returns function handles of the form u=@(x,y,z,t)
 
                          
