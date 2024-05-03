---
layout: default
title: linearWVIndexFromModeNumber
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 29
mathjax: true
---

#  linearWVIndexFromModeNumber

return the linear index into k_wv and l_wv from a mode number


---

## Declaration
```matlab
 index = linearWVIndexFromModeNumber(kMode,lMode,jMode)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer

## Returns
+ `linearIndex`  a non-negative integer number

## Discussion

  This function will return the linear index into the (k_wv,l_wv) arrays,
  given the mode numbers (kMode,lMode). Note that this will
  *not* normalize the mode to the primary mode number, but will
  throw an error.
 
          
