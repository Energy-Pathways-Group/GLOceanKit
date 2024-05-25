---
layout: default
title: setGeostrophicModes
parent: WVTransformBoussinesq
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  setGeostrophicModes

set amplitudes of the given geostrophic modes


---

## Declaration
```matlab
 [k,l] = setGeostrophicModes(self)
```
## Parameters
+ `kMode`  (optional) integer index, (k0 > -Nx/2 && k0 < Nx/2)
+ `lMode`  (optional) integer index, (l0 > -Ny/2 && l0 < Ny/2)
+ `jMode`  (optional) integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
+ `phi`  (optional) phase in radians, (0 <= phi <= 2*pi)
+ `u`  (optional) fluid velocity u (m/s)

## Returns
+ `k`  wavenumber k of the waves (radians/m)
+ `l`  wavenumber l of the waves (radians/m)

## Discussion

  Overwrite any existing amplitudes to any existing amplitudes
 
                    
