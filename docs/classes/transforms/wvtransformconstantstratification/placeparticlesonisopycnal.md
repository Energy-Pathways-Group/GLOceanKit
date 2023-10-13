---
layout: default
title: PlaceParticlesOnIsopycnal
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 15
mathjax: true
---

#  PlaceParticlesOnIsopycnal

MAS 1/10/18 - added intext ('int' or 'both') to give option of using int vs. int+ext fields for rho_prime


---

## Discussion
Also added rhoIsopycnal as output.
  Any floats with the same value of z will be moved to the same
  isopycnal.
 
  interpolation should probably be 'spline'.
  tolerance is in meters, 1e-8 should be fine.
