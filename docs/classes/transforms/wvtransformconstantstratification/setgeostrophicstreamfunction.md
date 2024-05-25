---
layout: default
title: setGeostrophicStreamfunction
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 29
mathjax: true
---

#  setGeostrophicStreamfunction

set a geostrophic streamfunction


---

## Declaration
```matlab
 setGeostrophicStreamfunction(psi)
```
## Parameters
+ `psi`  function handle that takes three arguments, psi(X,Y,Z)

## Discussion

  Clears A0 by setting a geostrophic streamfunction
 
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
 
      
