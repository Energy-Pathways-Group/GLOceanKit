---
layout: default
title: normalizeGeostrophicModeProperties
parent: WVGeostrophicSolutionGroup
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  normalizeGeostrophicModeProperties

returns properties of a geostrophic solution relative to the primary mode number


---

## Declaration
```matlab
 [kMode,lMode,jMode,A,phi] = normalizeGeostrophicModeProperties(self,kMode,lMode,jMode,A,phi)
```
## Parameters
+ `kMode`  integer
+ `lMode`  integer
+ `jMode`  non-negative integer
+ `A`  real-valued amplitude (m)
+ `phi`  real-valued phase (radians)

## Returns
+ `kMode`  integer
+ `lMode`  integer
+ `jMode`  non-negative integer
+ `A`  real-valued amplitude (m)
+ `phi`  real-valued phase (radians)

## Discussion

  This function will return the primary mode numbers (k,l,j),
  given the any valid mode numbers (k,l,j) and adjust the
  amplitude (A) and phase (phi), if necessary.
 
                        
