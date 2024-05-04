---
layout: default
title: linearIndexFromModeNumber
parent: WVOrthogonalSolutionGroup
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  linearIndexFromModeNumber

return the linear index from the primary mode number


---

## Declaration
```matlab
 index = linearIndexFromModeNumber(kMode,lMode,jMode)
```
## Parameters
+ `kMode`  non-negative integer
+ `lMode`  non-negative integer
+ `jMode`  non-negative integer

## Returns
+ `linearIndex`  a non-negative integer number

## Discussion

  This function will return the linear index into the A0 array,
  given the primary mode numbers (k,l,j). Note that this will
  *not* normalize the mode to the primary mode number, but will
  throw an error.
 
            
