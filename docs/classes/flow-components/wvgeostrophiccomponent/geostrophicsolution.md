---
layout: default
title: geostrophicSolution
parent: WVGeostrophicComponent
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  geostrophicSolution

return a real-valued analytical solution of the geostrophic mode


---

## Declaration
```matlab
 solution = geostrophicSolution(kMode,lMode,jMode,A,phi,options)
```
## Parameters
+ `kMode`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
+ `A`  amplitude in m.
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `shouldAssumeConstantN`  (optional) default 1
+ `amplitudeIsMaxU`  (optional) default 0

## Returns
+ `u`  fluid velocity, u = @(x,y,z,t)
+ `v`  fluid velocity, v = @(x,y,z,t)
+ `w`  fluid velocity, w = @(x,y,z,t)
+ `eta`  isopycnal displacement, eta = @(x,y,z,t)
+ `p`  pressure, p = @(x,y,z,t)

## Discussion

  Returns function handles of the form u=@(x,y,z,t)
 
                            
