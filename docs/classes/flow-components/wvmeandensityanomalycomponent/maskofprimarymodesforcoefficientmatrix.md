---
layout: default
title: maskOfPrimaryModesForCoefficientMatrix
parent: WVMeanDensityAnomalyComponent
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  maskOfPrimaryModesForCoefficientMatrix

returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.


---

## Declaration
```matlab
 mask = maskOfPrimaryModesForCoefficientMatrix(coefficientMatrix)
```
## Parameters
+ `coefficientMatrix`  a WVCoefficientMatrix type

## Returns
+ `mask`  matrix of size [Nk Nl Nj] with 1s and 0s

## Discussion

  Returns a 'mask' (matrix with 1s or 0s) indicating where
  different solution types live in the Ap, Am, A0 matrices.
 
        
Help for WVMeanDensityAnomalyComponent/maskOfPrimaryModesForCoefficientMatrix is inherited from superclass WVPrimaryFlowComponent
