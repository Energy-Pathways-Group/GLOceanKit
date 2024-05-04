---
layout: default
title: wvIndexFromModeNumber
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 44
mathjax: true
---

#  wvIndexFromModeNumber

return the linear index into k_wv and l_wv from a mode number


---

## Declaration
```matlab
 index = wvIndexFromModeNumber(kMode,lMode,jMode)
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
 
          
