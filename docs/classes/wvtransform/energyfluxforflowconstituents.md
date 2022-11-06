---
layout: default
title: energyFluxForFlowConstituents
parent: WVTransform
grand_parent: Classes
nav_order: 81
mathjax: true
---

#  energyFluxForFlowConstituents

returns the energy flux into each coefficient, from specific flow constituents


---

## Declaration
```matlab
 [Ep,Em,E0] = energyFluxForFlowConstituents(Uconstituent,gradUconstituent,options)
```
## Parameters
+ `Uconstituent`  `WVFlowConstituent` type for $$\vec{u} \cdot \nabla \vec{u}$$
+ `gradUconstituent`  `WVFlowConstituent` type for $$\vec{u} \cdot \nabla \vec{u}$$
+ `deltaT`  (optional) include the deltaT term in the Euler time step

## Returns
+ `Ep`  energy flux into the Ap coefficients
+ `Em`  energy flux into the Am coefficients
+ `E0`  energy flux into the A0 coefficients

## Discussion

                
