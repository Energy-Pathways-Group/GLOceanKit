---
layout: default
title: checkHermitian
parent: WVTransform
grand_parent: Classes
nav_order: 66
mathjax: true
---

#  checkHermitian

Check if the matrix is Hermitian. Report errors.


---

## Declaration
```matlab
 A = checkHermitian( A )
```
## Returns
+ `A`  matrix the same size as the input matrix

## Discussion

  This function makes assumptions about the structure of the matrix.
 
  The approach taken here is that the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
  numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
  conjugates. Also, the negative k wavenumbers for l=0.
 
      
