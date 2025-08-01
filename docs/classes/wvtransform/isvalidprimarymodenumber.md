---
layout: default
title: isValidPrimaryModeNumber
parent: WVTransform
grand_parent: Classes
nav_order: 108
mathjax: true
---

#  isValidPrimaryModeNumber

returns a boolean indicating whether (k,l,j) is a valid primary (non-conjugate) mode number


---

## Declaration
```matlab
 index = isValidPrimaryModeNumber(kMode,lMode,jMode)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer
+ `jMode`  non-negative integer

## Returns
+ `index`  a non-negative integer

## Discussion

  returns a boolean indicating whether (k,l,j) is a valid
  non-conjugate mode number according to how the property
  conjugateDimension is set.
 
            
