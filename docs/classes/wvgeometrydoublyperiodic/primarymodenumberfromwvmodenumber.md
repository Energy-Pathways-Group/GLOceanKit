---
layout: default
title: primaryModeNumberFromWVModeNumber
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 35
mathjax: true
---

#  primaryModeNumberFromWVModeNumber

takes any valid WV mode number and returns the primary mode number


---

## Declaration
```matlab
 [kMode,lMode] = primaryModeNumberFromWVModeNumber(kMode,lMode)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer

## Returns
+ `kMode`  integer
+ `lMode`  integer

## Discussion

  The function first confirms that the mode numbers are valid,
  and then converts any conjugate mode numbers to primary mode
  numbers.
 
  The result is affected by the chosen conjugateDimension.
 
            
