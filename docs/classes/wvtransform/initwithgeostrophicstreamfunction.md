---
layout: default
title: initWithGeostrophicStreamfunction
parent: WVTransform
grand_parent: Classes
nav_order: 112
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
 
  Clears variables Ap,Am,A0 and then sets the geostrophic streamfunction
      
