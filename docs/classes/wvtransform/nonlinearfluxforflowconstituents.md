---
layout: default
title: nonlinearFluxForFlowConstituents
parent: WVTransform
grand_parent: Classes
nav_order: 136
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

  This computes the nonlinear flux that results from a subset of flow
  constituents. The masks are applied to the coefficients Ap,Am,A0 before
  computing the nonlinear flux, $$\vec{u} \cdot \nabla \vec{u}$$. This
  function calls -nonlinearFluxWithGradientMasks.
 
              
