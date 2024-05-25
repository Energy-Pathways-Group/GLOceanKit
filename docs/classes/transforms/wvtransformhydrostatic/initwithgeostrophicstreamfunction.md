---
layout: default
title: initWithGeostrophicStreamfunction
parent: WVTransformHydrostatic
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  initWithGeostrophicStreamfunction

initialize with a geostrophic streamfunction


---

## Declaration
```matlab
 initWithGeostrophicStreamfunction(psi)
```
## Parameters
+ `psi`  function handle that takes three arguments, psi(X,Y,Z)

## Discussion

  Clears variables Ap,Am,A0 and then sets the geostrophic
  streamfunction.
 
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
 
      
