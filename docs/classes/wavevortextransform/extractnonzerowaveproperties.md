---
layout: default
title: extractNonzeroWaveProperties
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 85
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

## Discussion

  This function makes assumptions about the structure of the matrix.
        - Returns A: amplitude
  - Returns phi: phase
  - Returns linearIndex: linear index of matrix component
