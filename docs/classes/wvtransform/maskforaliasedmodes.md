---
layout: default
title: maskForAliasedModes
parent: WVTransform
grand_parent: Classes
nav_order: 122
mathjax: true
---

#  maskForAliasedModes

returns a mask with locations of modes that will alias with a quadratic multiplication.


---

## Declaration
```matlab
 AntiAliasFilter = maskForAliasedModes(self,options)
```
## Parameters
+ `jFraction`  (optional) fraction of vertical mode to assume are not aliased (default 2/3)

## Returns
+ `AntiAliasFilter`  mask aliased mode

## Discussion

  Returns a 'mask' (matrices with 1s or 0s) indicating where aliased wave
  modes are, assuming the 2/3 anti-aliasing rule for quadratic
  interactions. The reality is that the vertical modes do not follow this
  rule.
 
  Basic usage,
  AntiAliasMask = wvm.maskForAliasedModes();
  will return a mask that contains 1 at the locations of modes that will
  alias with a quadratic multiplication.
 
        
