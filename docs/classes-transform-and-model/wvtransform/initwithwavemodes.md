---
layout: default
title: initWithWaveModes
parent: WVTransform
grand_parent: Classes
nav_order: 109
mathjax: true
---

#  initWithWaveModes

initialize with the given wave modes


---

## Declaration
```matlab
 [omega,k,l] = initWithWaveModes( kMode, lMode, jMode, phi, Amp, signs)
```
## Parameters
+ `k`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `l`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `j`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `u`  fluid velocity u (m/s)
+ `sign`  sign of the frequency, +1 or -1

## Returns
+ `omega`  frequencies of the waves (radians/s)
+ `k`  wavenumber k of the waves (radians/m)
+ `l`  wavenumber l of the waves (radians/m)

## Discussion

  $$
  sin(k*x+l*y)*F_j*sin(omega*t + phi)
  $$
 
  Clears variables Ap,Am,A0 and then sets the given wave modes.
                      
