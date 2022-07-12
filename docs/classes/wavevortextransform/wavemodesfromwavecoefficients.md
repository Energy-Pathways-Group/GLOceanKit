---
layout: default
title: waveModesFromWaveCoefficients
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 178
mathjax: true
---

#  waveModesFromWaveCoefficients

This returns the properties of the waves being used in the


---

## Discussion
gridded simulation, as their properly normalized individual
  wave components. Very useful for debugging.
 
  Note that A_plus and A_minus each have half the inertial
  energy. This can be misleading, but the phasing is chosen to
  make it work. Never-the-less, we double/zero that component.
