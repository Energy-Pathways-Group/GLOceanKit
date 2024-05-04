---
layout: default
title: maskForAliasedModes
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 31
mathjax: true
---

#  maskForAliasedModes

returns a mask with locations of modes that will alias with a quadratic multiplication.


---

## Declaration
```matlab
 antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(Nx,Ny,Nz);
```
## Parameters
+ `Nx`  grid points in the x-direction
+ `Ny`  grid points in the y-direction
+ `Nz`  grid points in the z-direction (defuault 1)

## Returns
+ `antialiasMask`  mask aliased mode

## Discussion

  Returns a 'mask' (matrices with 1s or 0s) indicating where aliased wave
  modes are, assuming the 2/3 anti-aliasing rule for quadratic
  interactions.
 
  Technically one needs only restrict to 2/3s in each
  wavenumber direction. However, we prefer to maintain an
  isotropic effective grid size and instead restrict to a
  circle.
 
  Basic usage,
  ```matlab
  antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(8,8);
  ```
  will return a mask that contains 1 at the locations of modes that will
  alias with a quadratic multiplication.
 
            
