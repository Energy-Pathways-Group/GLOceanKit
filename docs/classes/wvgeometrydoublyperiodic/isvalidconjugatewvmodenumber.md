---
layout: default
title: isValidConjugateWVModeNumber
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  isValidConjugateWVModeNumber

return a boolean indicating whether (k,l) is a valid conjugate WV mode number


---

## Declaration
```matlab
 bool = isValidConjugateWVModeNumber(kMode,lMode)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer

## Returns
+ `bool`  [0 1]

## Discussion

  returns a boolean indicating whether (k,l) is a valid
  *conjugate* WV mode number. Even if a mode number is
  available in the DFT matrix, it does not mean it is a valid
  WV mode number, e.g., it may be removed due to aliasing.
 
  The result is affected by the chosen conjugateDimension.
 
  Any valid self-conjugate modes (i.e., k=l=0) will return 1.
 
          
