---
layout: default
title: transformFromDFTGridToWVGrid
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 41
mathjax: true
---

#  transformFromDFTGridToWVGrid

convert from DFT to WV grid


---

## Declaration
```matlab
 Azkl = transformFromDFTGridToWVGrid(self,Aklz)
```
## Parameters
+ `Aklz`  DFT format matrix of size [Nk_dft Nl_dft Nz] (equivalently [Nx Ny Nz]) where Nz can be of any length

## Returns
+ `Azkl`  WV format matrix of size [Nz Nkl_wv]

## Discussion

  This function will reformat the memory layout of a data
  structure on a DFT grid to one on a WV grid. The WV grid will
  respect the conditions set when this class was initialized
  (shouldAntialias, shouldExcludeNyquist,
  shouldExcludeConjugates).
 
  This function is not the fastest way to reformat your data.
  If high performance is required, you should 
 
        
