---
layout: default
title: normalizeWaveModeProperties
parent: WVInternalGravityWaveSolutionGroup
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  normalizeWaveModeProperties

returns properties of a internal gravity wave solutions relative to the primary mode number


---

## Declaration
```matlab
 [kMode,lMode,jMode,A,phi] = normalizeGeostrophicModeProperties(self,kMode,lMode,jMode,A,phi)
```
## Parameters
+ `kMode`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes)
+ `A`  amplitude in m/s.
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `omegasign`  sign of omega, [-1 1]

## Returns
+ `kMode`  integer index
+ `lMode`  integer index
+ `jMode`  integer index
+ `A`  amplitude in m.
+ `phi`  phase in radians
+ `omegasign`  sign of omega, [-1 1]

## Discussion

  This function will return the primary mode numbers (k,l,j),
  given the any valid mode numbers (k,l,j) and adjust the
  amplitude (A) and phase (phi), if necessary.
 
                            
