---
layout: default
title: addGeostrophicModes
parent: WVTransform
grand_parent: Classes
nav_order: 61
mathjax: true
---

#  addGeostrophicModes

add amplitudes of the given geostrophic modes


---

## Declaration
```matlab
 [k,l] = addGeostrophicModes(self)
```
## Parameters
+ `k`  integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `l`  integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `j`  integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
+ `phi`  phase in radians, (0 <= phi <= 2*pi)
+ `u`  fluid velocity u (m/s)

## Returns
+ `k`  wavenumber k of the waves (radians/m)
+ `l`  wavenumber l of the waves (radians/m)

## Discussion

  Add new amplitudes to any existing amplitudes
                  
