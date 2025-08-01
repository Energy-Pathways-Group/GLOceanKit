---
layout: default
title: enstrophyFactorForA0
parent: WVPrimaryFlowComponent
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  enstrophyFactorForA0

returns the enstrophy multiplier for the A0 coefficient matrix.


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
  produce the total enstrophy.
 
        
