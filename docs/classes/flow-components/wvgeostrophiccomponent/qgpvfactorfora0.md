---
layout: default
title: qgpvFactorForA0
parent: WVGeostrophicComponent
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  qgpvFactorForA0

returns the qgpv multiplier for the coefficient matrix.


---

## Declaration
```matlab
 totalEnergyFactor = totalEnergyFactorForCoefficientMatrix(coefficientMatrix)
```
## Parameters
+ `coefficientMatrix`  a WVCoefficientMatrix type

## Returns
+ `mask`  matrix of size [Nj Nkl]

## Discussion

  Returns a matrix of size wvt.spectralMatrixSize that
  multiplies the squared absolute value of this matrix to
  produce the total energy.
 
        
