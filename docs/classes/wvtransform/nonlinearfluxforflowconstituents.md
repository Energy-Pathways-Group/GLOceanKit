---
layout: default
title: nonlinearFluxForFlowConstituents
parent: WVTransform
grand_parent: Classes
nav_order: 125
mathjax: true
---

#  nonlinearFluxForFlowConstituents

returns the flux of each coefficient as determined by the nonlinear flux operation


---

## Declaration
```matlab
 [Fp,Fm,F0] = nonlinearFluxForFlowConstituents(Uconstituent,gradUconstituent)
```
## Parameters
+ `Uconstituent`  `WVFlowConstituent` type for $$\vec{u} \cdot \nabla \vec{u}$$
+ `gradUconstituent`  `WVFlowConstituent` type for $$\vec{u} \cdot \nabla \vec{u}$$

## Returns
+ `Fp`  flux into the Ap coefficients
+ `Fm`  flux into the Am coefficients
+ `F0`  flux into the A0 coefficients

## Discussion

              
