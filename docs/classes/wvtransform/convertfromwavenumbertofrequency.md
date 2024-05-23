---
layout: default
title: convertFromWavenumberToFrequency
parent: WVTransform
grand_parent: Classes
nav_order: 74
mathjax: true
---

#  convertFromWavenumberToFrequency

Summary


---

## Declaration
```matlab
 [varargout] = wvt.convertFromWavenumberToFrequency
```
## Parameters
+ `varargin`  WVT

## Returns
+ `varargout`  energyFrequency has dimensions $$(j,omegaVector)$$

## Discussion
This method transforms an instantaneous WVT field that is in the 
   horizontal wave number domain to the frequency domain.
 
   At this point the method is only set to return the total wave energy.
   In the future I plan to include an option for the user
   to transform any WVT field (Leticia)
 
        
