---
layout: default
title: compute
parent: WVOperation
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  compute

the promised variable


---

## Declaration
```matlab
 varargout = compute(wvt,varargin)
```
## Parameters
+ `wvt`  A WVTransform instance from which to compute the variable

## Returns
+ `varargout`  cell array of returned variables

## Discussion

  The compute operation is given the current state of the ocean
  (wvt) during each call, and it is expected to compute the
  variables from that state. You can, of course, use other 
  instance variables from your own custom subclass to
  implement computations with other dependencies.
  
        
