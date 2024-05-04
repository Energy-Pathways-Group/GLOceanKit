---
layout: default
title: maskForPrimaryCoefficients
parent: WVMeanDensityAnomalySolutionGroup
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  maskForPrimaryCoefficients

returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.


---

## Declaration
```matlab
 mask = maskForPrimaryCoefficients(coefficientMatrix)
```
## Parameters
+ `coefficientMatrix`  a WVCoefficientMatrix type

## Returns
+ `mask`  matrix of size [Nk Nl Nj] with 1s and 0s

## Discussion

  Returns a 'mask' (matrix with 1s or 0s) indicating where
  different solution types live in the Ap, Am, A0 matrices.
 
        
