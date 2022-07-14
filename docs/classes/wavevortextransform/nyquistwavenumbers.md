---
layout: default
title: nyquistWavenumbers
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 120
mathjax: true
---

#  nyquistWavenumbers

Returns a matrix with 1s at the Nyquist frequencies.


---

## Declaration
```matlab
 A = nyquistWavenumbers( A )
```
## Discussion

  This function makes assumptions about the structure of the matrix.
 
  Returns a matrix of the same size as A with 1s at the Nyquist
  frequencies.
 
      - Returns A: matrix the same size as the input matrix
