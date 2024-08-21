---
layout: default
title: effectiveHorizontalGridResolution
parent: WVTransform
grand_parent: Classes
nav_order: 84
mathjax: true
---

#  effectiveHorizontalGridResolution

returns the effective grid resolution in meters


---

## Declaration
```matlab
 flag = effectiveHorizontalGridResolution(other)
```
## Returns
+ `effectiveHorizontalGridResolution`  double

## Discussion

  The effective grid resolution is the highest fully resolved
  wavelength in the model. This value takes into account
  anti-aliasing, and is thus appropriate for setting damping
  operators.
 
      
