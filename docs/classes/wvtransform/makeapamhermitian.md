---
layout: default
title: makeApAmHermitian
parent: WVTransform
grand_parent: Classes
nav_order: 152
mathjax: true
---

#  makeApAmHermitian

Forces the Ap/Am matrices to have the correct symmetries


---

## Declaration
```matlab
 [Ap,Am] = makeApAmHermitian(Ap,Am)
```
## Returns
+ `Ap`  matrix the same size as the input matrix
+ `Am`  matrix the same size as the input matrix

## Discussion

  This function is NOT a true "Make Hermitian" function because it
  the Ap/Am matrices do not require k=l=0 to be real.
 
  If conjugateDimension == 2, then the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
  numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
  conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
  numbers are set to zero to avoid complications.
 
        
