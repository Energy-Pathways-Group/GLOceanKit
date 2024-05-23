---
layout: default
title: enstrophyFactorForA0
parent: WVGeostrophicComponent
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  enstrophyFactorForA0

returns the qgpv multiplier for the A0 coefficient matrix.


---

## Declaration
```matlab
 qgpvFactor = qgpvFactorForA0()
```
## Returns
+ `qgpvFactor`  matrix of size [Nj Nkl]

## Discussion

  Returns a matrix of size wvt.spectralMatrixSize that
  multiplies the A0 matrix so that when transformed with the Fg
  modes will return QGPV.
 
      
