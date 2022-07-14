---
layout: default
title: maskForAliasedModes
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 119
mathjax: true
---

#  maskForAliasedModes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


---

## Discussion

  Returns a 'mask' (matrices with 1s or 0s) indicating where aliased wave
  modes are, assuming the 2/3 anti-aliasing rule for quadratic
  interactions. The reality is that the vertical modes do not follow this
  rule.
 
  Basic usage,
  AntiAliasMask = wvm.maskForAliasedModes();
  will return a mask that contains 1 at the locations of modes that will
  alias with a quadratic multiplication.
