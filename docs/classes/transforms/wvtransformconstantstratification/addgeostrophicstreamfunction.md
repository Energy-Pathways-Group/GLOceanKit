---
layout: default
title: addGeostrophicStreamfunction
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 16
mathjax: true
---

#  addGeostrophicStreamfunction

add a geostrophic streamfunction to existing geostrophic motions


---

## Declaration
```matlab
 addGeostrophicStreamfunction(psi)
```
## Parameters
+ `psi`  function handle that takes three arguments, psi(X,Y,Z)

## Discussion

  The geostrophic streamfunction is added to the existing values in `A0`
 
  The geostrophic streamfunction, $$\psi$$, is defined such that
 
  $$
  u= - \frac{\partial \psi}{\partial y}
  $$
 
  $$
  v=\frac{\partial \psi}{\partial x}
  $$
 
  $$
  N^2 \eta = \frac{g}{\rho_0} \rho = - f \frac{\partial \psi}{\partial z}
  $$
 
  Note that a streamfunction also projects onto the
  mean-density-anomaly (MDA) component of the flow, and thus it
  is not strictly geostrophic.
 
      
