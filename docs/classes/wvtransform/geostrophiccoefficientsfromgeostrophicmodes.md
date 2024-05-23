---
layout: default
title: geostrophicCoefficientsFromGeostrophicModes
parent: WVTransform
grand_parent: Classes
nav_order: 99
mathjax: true
---

#  geostrophicCoefficientsFromGeostrophicModes

Returns the indices (and re-normalized values) of the geostropic mode appropriate for the A0 matrix.


---

## Declaration
```matlab
 [kIndex,lIndex,jIndex,A0Amp] = geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, phi, u, signs)
```
## Parameters
+ `kMode`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `u`  fluid velocity u (m/s)

## Discussion

  Returns the indices (and re-normalized values) of the geostrophic mode
  appropriate for the A0 matrices. This works in conjunction with the
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
  u is the fluid velocity U
 
              
