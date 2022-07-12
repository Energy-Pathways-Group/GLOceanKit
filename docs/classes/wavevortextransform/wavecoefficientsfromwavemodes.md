---
layout: default
title: waveCoefficientsFromWaveModes
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 176
mathjax: true
---

#  waveCoefficientsFromWaveModes

Returns the indices (and re-normalized values) of the wave mode


---

## Discussion
appropriate for the Ap, Am matrices. This works in conjunction with the
  makeHermitian function, which then sets the appropriate conjugate. At the
  moment we made the (perhaps bad) choice that the negative l components
  are redundant, but to take advantage of the FFT, we may change this in
  the future.
  
  For example, wave mode with l<0, is equivalent to a wave mode with l>0
  and the signs fliipped on all the other quantities.
 
  The values given must meet the following requirements:
  (k0 > -Nx/2 && k0 < Nx/2)
  (l0 > -Ny/2 && l0 < Ny/2)
  (j0 >= 1 && j0 <= nModes)
  phi is in radians, from 0-2pi
  Amp is the fluid velocity U
  sign is +/-1, indicating the sign of the frequency.
