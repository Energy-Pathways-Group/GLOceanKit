---
layout: default
title: addGeostrophicStreamfunction
parent: WVTransform
grand_parent: Classes
nav_order: 60
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

  The geostrophic streamfunction, $\psi$, is defined such that
  $$
  u= - \frac{\partial \psi}{\partial y}
  $$
  
  $$
  v=\frac{\partial \psi}{\partial x}
  $$
  
  $$
  N^2 \eta = - f \frac{\partial \psi}{\partial z}
  $$
 
  The geostrophic streamfunction is added to the existing values in `A0`
      
