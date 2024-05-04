---
layout: default
title: subscriptIndicesFromPrimaryModeNumber
parent: WVOrthogonalSolutionGroup
grand_parent: Classes
nav_order: 15
mathjax: true
---

#  subscriptIndicesFromPrimaryModeNumber

return subscript indices for a given mode number


---

## Declaration
```matlab
 [kIndex,lIndex,jIndex] = subscriptIndicesFromPrimaryModeNumber(kMode,lMode,jMode)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer
+ `jMode`  non-negative integer

## Returns
+ `kIndex`  a positive integer
+ `lIndex`  a positive integer
+ `jIndex`  a positive integer

## Discussion

  This function will return the subscript indices into any
  coefficient matrix, given the mode numbers (k,l,j). Note that
  this will *not* normalize the mode to the primary mode
  number, but will throw an error if you request the conjugate
  mode number.
 
                
