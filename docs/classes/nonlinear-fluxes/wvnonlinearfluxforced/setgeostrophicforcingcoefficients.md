---
layout: default
title: setGeostrophicForcingCoefficients
parent: WVNonlinearFluxForced
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  setGeostrophicForcingCoefficients

set forcing values for the geostrophic part of the flow


---

## Declaration
```matlab
 setGeostrophicForcingCoefficients(A0bar,options)
```
## Parameters
+ `A0bar`  A0 'mean' value to relax to
+ `MA0`  (optional) forcing mask, A0. 1s at the forced modes, 0s at the unforced modes. If it is left blank, then it will be produced using the nonzero values of A0bar
+ `tau0`  (optional) relaxation time

## Discussion

  Forcing takes the following form,
 
  $$
  \frac{\partial}{\partial t} A_0^{klj} = \underbrace{M_{A_0}^{klj} \left(\bar{A}_0^{klj}  - A_0^{klj} \right)/ \tau_0}_{F_\textrm{force}} + F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
  $$
 
  where $$M_{A_0}^{klj}$$ are masks (1s and 0s),
  $$\bar{A}_0^{klj}$$ are mean amplitudes, and $$\tau_0$$
  are time scales. If the time scale is set to 0, then the mean
  amplitudes remain fixed for all time. In that limit, the
  equations can be written as,
 
  $$
  \frac{\partial}{\partial t} A_0^{klj} = \neg M_{A_0}^{klj} \left( F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj} \right)
  $$
 
          
