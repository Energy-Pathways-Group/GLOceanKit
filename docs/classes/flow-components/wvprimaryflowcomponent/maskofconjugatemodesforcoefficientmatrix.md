---
layout: default
title: maskOfConjugateModesForCoefficientMatrix
parent: WVPrimaryFlowComponent
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  maskOfConjugateModesForCoefficientMatrix

returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.


---

## Declaration
```matlab
 mask = maskOfConjugateModesForCoefficientMatrix(self,coefficientMatrix)
```
## Parameters
+ `coefficientMatrix`  a WVCoefficientMatrix type

## Returns
+ `mask`  matrix of size [Nk Nl Nj] with 1s and 0s

## Discussion

  Returns a 'mask' (matrix with 1s or 0s) indicating where
  different solution types live in the Ap, Am, A0 matrices.
 
        
