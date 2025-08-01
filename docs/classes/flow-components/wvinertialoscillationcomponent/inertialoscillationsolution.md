---
layout: default
title: inertialOscillationSolution
parent: WVInertialOscillationComponent
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  inertialOscillationSolution

return a real-valued analytical solution of the internal gravity wave mode


---

## Declaration
```matlab
 solution = internalGravityWaveSolution(self,kMode,lMode,jMode,A,phi,omegasign,options)
```
## Parameters
+ `kMode`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
+ `A`  amplitude in m/s.
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `shouldAssumeConstantN`  (optional) default 1

## Returns
+ `u`  fluid velocity, u = @(x,y,z,t)
+ `v`  fluid velocity, v = @(x,y,z,t)
+ `w`  fluid velocity, w = @(x,y,z,t)
+ `eta`  isopycnal displacement, eta = @(x,y,z,t)
+ `p`  pressure, p = @(x,y,z,t)

## Discussion

  Returns function handles of the form u=@(x,y,z,t)
 
                          
