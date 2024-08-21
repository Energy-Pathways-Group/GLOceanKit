---
layout: default
title: effectiveVerticalGridResolution
parent: WVTransformBoussinesq
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  effectiveVerticalGridResolution

returns the effective vertical grid resolution in meters


---

## Declaration
```matlab
 flag = effectiveVerticalGridResolution(other)
```
## Returns
+ `effectiveVerticalGridResolution`  double

## Discussion

  The effective grid resolution is the highest fully resolved
  wavelength in the model. This value takes into account
  anti-aliasing, and is thus appropriate for setting damping
  operators.
 
      
