---
layout: default
title: totalEnergyFactorForCoefficientMatrix
parent: WVInertialOscillationComponent
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  totalEnergyFactorForCoefficientMatrix

returns the total energy multiplier for the coefficient matrix.


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
 
        
