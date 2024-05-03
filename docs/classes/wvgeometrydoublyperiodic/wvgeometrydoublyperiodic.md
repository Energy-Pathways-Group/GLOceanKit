---
layout: default
title: WVGeometryDoublyPeriodic
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  WVGeometryDoublyPeriodic

create a geometry for a  doubly periodic domain


---

## Declaration
```matlab
  self = WVGeometryDoublyPeriodic(Lxy, Nxy, options)
```
## Parameters
+ `Lxy`  length of the domain (in meters) in the two periodic coordinate directions, e.g. [Lx Ly]
+ `Nxy`  number of grid points in the two coordinate directions, e.g. [Nx Ny]
+ `conjugateDimension`  (optional) set which dimension in the DFT grid is assumed to have the redundant conjugates (1 or 2), default is 2
+ `shouldAntialias`  (optional) set whether the WV grid excludes the quadratically aliased modes [0 1] (default 1)
+ `shouldExcludeNyquist`  (optional) set whether the WV grid excludes Nyquist modes[0 1] (default 1)
+ `shouldExludeConjugates`  (optional) set whether the WV grid excludes conjugate modes [0 1] (default 1)

## Returns
+ `geom`  a new WVGeometryDoublyPeriodic instance

## Discussion

                  
