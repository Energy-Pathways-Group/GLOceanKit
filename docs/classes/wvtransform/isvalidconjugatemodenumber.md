---
layout: default
title: isValidConjugateModeNumber
parent: WVTransform
grand_parent: Classes
nav_order: 104
mathjax: true
---

#  isValidConjugateModeNumber

returns a boolean indicating whether (k,l,j) is a valid conjugate mode number


---

## Declaration
```matlab
 index = isValidConjugateModeNumber(kMode,lMode,jMode)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer
+ `jMode`  non-negative integer

## Returns
+ `index`  a non-negative integer

## Discussion

  returns a boolean indicating whether (k,l,j) is a valid
  conjugate mode number according to how the property
  conjugateDimension is set.
 
  Any valid self-conjugate modes (i.e., k=l=0) will return 1.
 
            
