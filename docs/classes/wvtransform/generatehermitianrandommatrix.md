---
layout: default
title: generateHermitianRandomMatrix
parent: WVTransform
grand_parent: Classes
nav_order: 90
mathjax: true
---

#  generateHermitianRandomMatrix

Generate a 3D matrix to be Hermitian, except at k=l=0


---

## Declaration
```matlab
 A = generateHermitianRandomMatrix( size, options )
```
## Parameters
+ `size`  size of the matrix to generate
+ `shouldExcludeNyquist`  optional (default 1) will set Nyquist frequencies to zero if true.
+ `allowMeanPhase`  optional (default 0) will all a compex component to values at index (1,1,:) if set to true

## Returns
+ `A`  Hermitian conjugate matrix of given size

## Discussion

  This function makes assumptions about the structure of the matrix.
 
  Normally values at index (1,1,:) must be self-conjugate (and therefore
  real), but we added a flag to allow a phase, which is helpful for
  randomizing the phase of the inertial oscillations.
 
            
