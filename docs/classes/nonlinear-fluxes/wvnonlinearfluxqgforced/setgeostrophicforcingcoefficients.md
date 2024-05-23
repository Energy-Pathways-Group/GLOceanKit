---
layout: default
title: setGeostrophicForcingCoefficients
parent: WVNonlinearFluxQGForced
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  setGeostrophicForcingCoefficients

set forcing values for the geostrophic part of the flow


---

## Declaration
```matlab
 varargout = compute(wvt,varargin)
```
## Parameters
+ `A0bar`  A0 'mean' value to relax to
+ `MA0`  (optional) forcing mask, A0. 1s at the forced modes, 0s at the unforced modes. If it is left blank, then it will be produced using the nonzero values of A0bar
+ `tau0`  (optional) relaxation time

## Returns
+ `varargout`  cell array of returned variables

## Discussion

  $$
  \frac{\partial}{\partial t} A_0^{klj} = \underbrace{M_{A_0}^{klj} \left(\bar{A}_0^{klj}  - A_0^{klj} \right)/ \tau_0}_{F_\textrm{force}} + F_0^{klj} + F_\textrm{damp}^{klj}
  $$
 
            
