---
layout: default
title: waveModesFromWaveCoefficients
parent: WVTransform
grand_parent: Classes
nav_order: 216
mathjax: true
---

#  waveModesFromWaveCoefficients

Returns normalized amplitudes and phases of all waves


---

## Declaration
```matlab
 [omega, alpha, k, l, mode, phi, A, norm] = waveModesFromWaveCoefficients()
```
## Discussion

  This returns the properties of the waves being used in the
  gridded simulation, as their properly normalized individual
  wave components. Very useful for debugging.
 
  Note that A_plus and A_minus each have half the inertial
  energy. This can be misleading, but the phasing is chosen to
  make it work. Never-the-less, we double/zero that component.
 
    
