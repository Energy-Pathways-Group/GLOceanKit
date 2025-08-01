---
layout: default
title: transformToSpatialDomainAtPosition
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 45
mathjax: true
---

#  transformToSpatialDomainAtPosition

transform from $$(k,l)$$ on the DFT grid to $$(x,y)$$ at any position


---

## Declaration
```matlab
 u = transformToSpatialDomainAtPosition(u_bar)
```
## Parameters
+ `u_bar`  a complex-valued matrix of size [Nk_dft Nl_dft]

## Returns
+ `u`  a real-valued matrix of size [length(x) length(y)]

## Discussion

  Performs an inverse non-uniform Fourier transform to take a
  matrix from the DFT grid back to the spatial domain at any
  set of points (x,y).
 
        
