---
layout: default
title: setGeostrophicStreamfunction
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 33
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
 
  Consider a shallow eddy where the density anomaly sits close to the
  surface. This example was constructed in [Early, Hernández-Dueñas, Smith,
  and Lelong (2024)](https://arxiv.org/abs/2403.20269)
 
  ```matlab
  x0 = 3*Lx/4;
  y0 = Ly/2;
 
  Le = 80e3;
  He = 300;
  U = 0.20; % m/s
 
  H = @(z) exp(-(z/He/sqrt(2)).^2 );
  F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
  psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*H(z).*(F(x,y) - (pi*Le*Le/(wvt.Lx*wvt.Ly)));
 
  wvt.setGeostrophicStreamfunction(psi);
  ```
 
        
