---
layout: default
title: redundantHermitianCoefficients
parent: WVTransform
grand_parent: Classes
nav_order: 140
mathjax: true
---

#  redundantHermitianCoefficients

Returns a matrix with 1s at the 'redundant' hermiation indices.


---

## Declaration
```matlab
 A = redundantHermitianCoefficients( A )
```
## Returns
+ `A`  matrix the same size as the input matrix

## Discussion

  This function makes assumptions about the structure of the matrix.
 
  Returns a matrix the same size as A with 1s at the 'redundant'
  hermiation indices.
 
      
