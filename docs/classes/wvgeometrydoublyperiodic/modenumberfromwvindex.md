---
layout: default
title: modeNumberFromWVIndex
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 33
mathjax: true
---

#  modeNumberFromWVIndex

return mode number from a linear index into a WV matrix


---

## Declaration
```matlab
 [kMode,lMode] = modeNumberFromWVIndex(self,linearIndex)
```
## Parameters
+ `linearIndex`  a non-negative integer number

## Returns
+ `kMode`  integer
+ `lMode`  integer

## Discussion

  This function will return the mode numbers (kMode,lMode)
  given some linear index into a WV structured matrix.
 
          
