---
layout: default
title: transformFromSpatialDomain
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 41
mathjax: true
---

#  transformFromSpatialDomain

transform from $$(x,y,z)$$ to $$(k,l,z)$$ on the DFT grid


---

## Declaration
```matlab
 u_bar = transformFromSpatialDomain(u)
```
## Parameters
+ `u`  a real-valued matrix of size [Nx Ny Nz]

## Returns
+ `u_bar`  a complex-valued matrix of size [Nk_dft Nl_dft Nz]

## Discussion

  Performs a Fourier transform in the x and y direction. The
  resulting matrix is on the DFT grid.
 
        
