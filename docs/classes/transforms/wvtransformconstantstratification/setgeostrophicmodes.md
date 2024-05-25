---
layout: default
title: setGeostrophicModes
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 28
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
+ `kMode`  integer index, (kMode > -Nx/2 && kMode < Nx/2)
+ `lMode`  integer index, (lMode > -Ny/2 && lMode < Ny/2)
+ `j`  integer index, (j >= 1 && j <= nModes), unless k=l=j=0
+ `phi`  (optional) phase in radians, (0 <= phi <= 2*pi), default 0
+ `u`  fluid velocity u (m/s)

## Returns
+ `k`  wavenumber k of the kModes (radians/m)
+ `l`  wavenumber l of the lModes (radians/m)

## Discussion

  Set the amplitude of the given geostrophic modes by
  overwriting any existing amplitudes. The parameters are given
  as [horizontal and vertical modes](/users-guide/wavenumber-modes-and-indices.html),
  and the function will return the associated [horizontal wavenumbers](/users-guide/wavenumber-modes-and-indices.html)
  of those modes.
 
  For example,
 
  ```matlab
  wvt.addGeostrophicModes(kMode=0,lMode=1,jMode=1,u=0.5);
  ```
 
  will add a geostrophic mode.
 
                    
