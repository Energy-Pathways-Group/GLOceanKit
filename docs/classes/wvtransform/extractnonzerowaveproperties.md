---
layout: default
title: extractNonzeroWaveProperties
parent: WVTransform
grand_parent: Classes
nav_order: 92
mathjax: true
---

#  extractNonzeroWaveProperties

Takes a Hermitian matrix and returns the amplitude and phase of nonzero components


---

## Declaration
```matlab
 [A,phi,linearIndex] = ExtractNonzeroWaveProperties(Matrix)
```
## Parameters
+ `Matrix`  Hermitian conjugate matrix

## Returns
+ `A`  amplitude
+ `phi`  phase
+ `linearIndex`  linear index of matrix component

## Discussion

  This function makes assumptions about the structure of the matrix.
            
