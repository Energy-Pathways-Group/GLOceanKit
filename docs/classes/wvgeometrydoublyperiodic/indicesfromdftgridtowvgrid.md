---
layout: default
title: indicesFromDFTGridToWVGrid
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 16
mathjax: true
---

#  indicesFromDFTGridToWVGrid

indices to convert from DFT to WV grid


---

## Declaration
```matlab
 dftToWVIndices = indicesFromDFTGridToWVGrid(self,Nz)
```
## Parameters
+ `Nz`  length of the outer dimension (default 1)

## Returns
+ `dftToWVIndices`  indices into a DFT matrix

## Discussion

  This function returns indices to quickly reformat the memory
  layout of a data structure on a DFT grid to one on a WV grid.
  The resulting WV grid will respect the conditions set when
  this class was initialized (shouldAntialias,
  shouldExcludeNyquist, shouldExcludeConjugates).
 
  This function is should generally be faster than the function
  transformFromDFTGridToWVGrid if you cache these indices.
 
        
