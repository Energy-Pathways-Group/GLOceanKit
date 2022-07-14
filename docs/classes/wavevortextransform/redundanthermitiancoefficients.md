---
layout: default
title: redundantHermitianCoefficients
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 131
mathjax: true
---

#  redundantHermitianCoefficients

Returns a matrix with 1s at the 'redundant' hermiation indices.


---

## Declaration
```matlab
 A = redundantHermitianCoefficients( A )
```
## Discussion

  This function makes assumptions about the structure of the matrix.
 
  Returns a matrix the same size as A with 1s at the 'redundant'
  hermiation indices.
 
      - Returns A: matrix the same size as the input matrix
