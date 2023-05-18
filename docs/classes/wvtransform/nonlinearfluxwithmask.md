---
layout: default
title: nonlinearFluxWithMask
parent: WVTransform
grand_parent: Classes
nav_order: 128
mathjax: true
---

#  nonlinearFluxWithMask

returns the flux of each coefficient as determined by the nonlinear flux


---

## Declaration
```matlab
 [Fp,Fm,F0] = nonlinearFluxWithMask(mask)
```
## Parameters
+ `mask`  mask applied to all constituents

## Returns
+ `Fp`  flux into the Ap coefficients
+ `Fm`  flux into the Am coefficients
+ `F0`  flux into the A0 coefficients

## Discussion
operation and the given mask.
 
  The mask is applied to the coefficients Ap,Am,A0 before computing the
  nonlinear flux. This is useful for zeroing wavenumbers at given total
  wavenumber or frequency, for example.
 
  The nonlinear flux used is the unforced, invicid equations.
 
            
