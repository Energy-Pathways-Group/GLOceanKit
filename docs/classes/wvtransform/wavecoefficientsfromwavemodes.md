---
layout: default
title: waveCoefficientsFromWaveModes
parent: WVTransform
grand_parent: Classes
nav_order: 214
mathjax: true
---

#  waveCoefficientsFromWaveModes

Returns the indices (and re-normalized values) of the wave mode appropriate for the Ap, Am matrices.


---

## Declaration
```matlab
 [kIndex,lIndex,jIndex,ApAmp,AmAmp] = waveCoefficientsFromWaveModes(kMode, lMode, jMode, phi, u, signs)
```
## Parameters
+ `kMode`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `Amp`  fluid velocity u (m/s)
+ `sign`  sign of the frequency, +1 or -1

## Discussion

  Returns the indices (and re-normalized values) of the wave mode
  appropriate for the Ap, Am matrices. This works in conjunction with the
  makeHermitian function, which then sets the appropriate conjugate. At the
  moment we made the (perhaps bad) choice that the negative l components
  are redundant, but to take advantage of the FFT, we may change this in
  the future.
  
  For example, wave mode with l<0, is equivalent to a wave mode with l>0
  and the signs flipped on all the other quantities.
 
  The values given must meet the following requirements:
  (k0 > -Nx/2 && k0 < Nx/2)
  (l0 > -Ny/2 && l0 < Ny/2)
  (j0 >= 1 && j0 <= nModes)
  phi is in radians, from 0-2pi
  Amp is the fluid velocity U
  sign is +/-1, indicating the sign of the frequency.
 
                
