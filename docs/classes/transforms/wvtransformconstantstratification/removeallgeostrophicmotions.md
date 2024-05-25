---
layout: default
title: removeAllGeostrophicMotions
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 25
mathjax: true
---

#  removeAllGeostrophicMotions

remove all geostrophic motions


---

## Declaration
```matlab
 removeAllGeostrophicMotions()
```
## Discussion

  All geostrophic motions are removed by setting A0 to zero.
 
  **Note** that this does *not* remove the mean density anomaly
  (mda) part of the solution, just the geostrophic part. Thus,
  this function will not clear all parts of a geostrophic
  streamfunction.
 
    
