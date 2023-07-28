---
layout: default
title: energyFluxFromNonlinearFlux
parent: WVTransform
grand_parent: Classes
nav_order: 80
mathjax: true
---

#  energyFluxFromNonlinearFlux

converts nonlinear flux into energy flux


---

## Declaration
```matlab
 [Ep,Em,E0] = energyFluxFromNonlinearFlux(Fp,Fm,F0,options)
```
## Parameters
+ `Fp`  nonlinear flux into the Ap coefficients
+ `Fm`  nonlinear flux into the Am coefficients
+ `F0`  nonlinear flux into the A0 coefficients
+ `deltaT`  (optional) include the deltaT term in the Euler time step

## Returns
+ `Ep`  energy flux into the Ap coefficients
+ `Em`  energy flux into the Am coefficients
+ `E0`  energy flux into the A0 coefficients

## Discussion

  Multiplies the nonlinear flux (Fp,Fm,F0) by the appropriate coefficients
  to convert into an energy flux.
 
  Optional parameter deltaT added by Bailey Avila: This equation is C17 in
  the manuscript, but with addition of 1st term on LHS of C16 converted to
  energy using Apm_TE_factor or A0_TE_factor This differs from energyFlux
  due to the importance of the 2*F*F*deltaT in equation C16 at the initial
  condition.
 
                  
