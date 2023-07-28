---
layout: default
title: makeHermitian
parent: WVTransform
grand_parent: Classes
nav_order: 121
mathjax: true
---

#  makeHermitian

Forces a 3D matrix to be Hermitian


---

## Declaration
```matlab
 A = makeHermitian( A )
```
## Returns
+ `A`  matrix the same size as the input matrix

## Discussion

  This function makes assumptions about the structure of the matrix.
 
  The approach taken here is that the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
  numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
  conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
  numbers are set to zero to avoid complications.
 
  This function is NOT a true "Make Hermitian" function because it
  doesn't force the k=l=0 to be real.
 
      
