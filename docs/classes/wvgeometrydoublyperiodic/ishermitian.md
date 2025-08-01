---
layout: default
title: isHermitian
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 19
mathjax: true
---

#  isHermitian

Check if the matrix is Hermitian. Report errors.


---

## Declaration
```matlab
 isHermitian( A )
```
## Parameters
+ `A`  matrix of size [K L Z]
+ `shouldReportErrors`  flag indicating whether or not error should be reported

## Returns
+ `bool`  flag indicating whether or not the matrix is Hermitian

## Discussion

  This algorithm checks whether any 2 or 3 dimensional matrix
  is Hermitian in the first two dimensions, following the data
  structure of a 2D DFT algorithm. The third dimension can be
  any length, including length 1.
 
  Errors can be reported indicating which entries are not
  conjugate.
 
  This algorithm is not fast.
 
          
