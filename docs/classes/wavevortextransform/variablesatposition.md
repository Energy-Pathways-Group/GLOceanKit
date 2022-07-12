---
layout: default
title: variablesAtPosition
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 172
---

#  variablesAtPosition

Primary method for accessing the dynamical variables on the


---

## Discussion
at any position or time.
 
  The method argument specifies how off-grid values should be
  interpolated. Use 'exact' for the slow, but accurate,
  spectral interpolation. Otherwise use 'spline' or some other
  method used by Matlab's interp function.
 
  Valid variable options are 'u', 'v', 'w', 'rho_prime', and
  'zeta'.
