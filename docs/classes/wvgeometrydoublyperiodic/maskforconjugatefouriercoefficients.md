---
layout: default
title: maskForConjugateFourierCoefficients
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 31
mathjax: true
---

#  maskForConjugateFourierCoefficients

a mask indicate the components that are redundant conjugates


---

## Declaration
```matlab
 mask = WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(Nx,Ny,options);
```
## Parameters
+ `Nx`  grid points in the x-direction
+ `Ny`  grid points in the y-direction
+ `conjugateDimension`  (optional) set which dimension in the DFT grid is assumed to have the redundant conjugates (1 or 2), default is 2

## Returns
+ `matrix`  matrix containing linear indices

## Discussion

            
