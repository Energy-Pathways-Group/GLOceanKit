---
layout: default
title: transformFromWVGridToDFTGrid
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 43
mathjax: true
---

#  transformFromWVGridToDFTGrid

convert from a WV to DFT grid


---

## Declaration
```matlab
 Aklz = transformFromWVGridToDFTGrid(self,Azkl)
```
## Parameters
+ `Azkl`  WV format matrix of size [Nz Nkl_wv] where Nz can be of any length
+ `isHalfComplex`  (optional) set whether the DFT grid excludes modes iL>Ny/2 [0 1] (default 0)

## Returns
+ `Aklz`  DFT format matrix of size [Nk_dft Nl_dft Nz] (equivalently [Nx Ny Nz])

## Discussion

  This function will reformat the memory layout of a data
  structure on a WV grid to one on a DFT grid. If the option
  isHalfComplex is selected, then it will not set values for
  iL>Ny/2, which are ignored by a 'symmetric' fft.
 
          
