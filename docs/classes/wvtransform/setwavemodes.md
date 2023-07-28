---
layout: default
title: setWaveModes
parent: WVTransform
grand_parent: Classes
nav_order: 160
mathjax: true
---

#  setWaveModes

set amplitudes of the given wave modes


---

## Declaration
```matlab
 [omega,k,l] = setWaveModes(kMode, lMode, jMode, phi, u, signs)
```
## Parameters
+ `kMode`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `Amp`  fluid velocity u (m/s)
+ `sign`  sign of the frequency, +1 or -1

## Returns
+ `omega`  frequencies of the waves (radians/s)
+ `k`  wavenumber k of the waves (radians/m)
+ `l`  wavenumber l of the waves (radians/m)

## Discussion

  Overwrite any existing wave modes with the given new values
                      
