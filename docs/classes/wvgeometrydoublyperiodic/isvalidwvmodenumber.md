---
layout: default
title: isValidWVModeNumber
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  isValidWVModeNumber

return a boolean indicating whether (k,l) is a valid WV mode number


---

## Declaration
```matlab
 bool = isValidWVModeNumber(kMode,lMode)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer

## Returns
+ `bool`  [0 1]

## Discussion

  returns a boolean indicating whether (k,l) is a valid WV mode
  number. Even if a mode number is available in the DFT matrix,
  it does not mean it is a valid WV mode number, e.g., it may
  be removed due to aliasing.
 
          
