---
layout: default
title: normalizeGeostrophicModeProperties
parent: WVGeostrophicSolutionGroup
grand_parent: Classes
nav_order: 13
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
+ `kMode`  non-negative integer
+ `lMode`  non-negative integer
+ `jMode`  non-negative integer

## Returns
+ `kIndex`  a positive integer
+ `lIndex`  a positive integer
+ `jIndex`  a positive integer

## Discussion

  This function will return the primary mode numbers (k,l,j),
  given the any valid mode numbers (k,l,j) and adjust the
  amplitude (A) and phase (phi), if necessary.
 
                
