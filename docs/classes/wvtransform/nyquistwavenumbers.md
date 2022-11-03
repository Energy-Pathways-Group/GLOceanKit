---
layout: default
title: nyquistWavenumbers
parent: WVTransform
grand_parent: Classes
nav_order: 125
mathjax: true
---

#  nyquistWavenumbers

Returns a matrix with 1s at the Nyquist frequencies.


---

## Declaration
```matlab
 A = nyquistWavenumbers( A )
```
## Returns
+ `A`  matrix the same size as the input matrix

## Discussion

  This function makes assumptions about the structure of the matrix.
 
  Returns a matrix of the same size as A with 1s at the Nyquist
  frequencies.
 
      
