---
layout: default
title: addMeanDensityAnomaly
parent: WVTransformHydrostatic
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  addMeanDensityAnomaly

add inertial motions to existing inertial motions


---

## Declaration
```matlab
 addMeanDensityAnomaly(eta)
```
## Parameters
+ `eta`  function handle that takes a single argument, eta(Z)

## Discussion

  The amplitudes of the inertial part of the flow will be added
  to the existing inertial part of the flow.
 
  ```matlab
  U_io = 0.2;
  Ld = wvt.Lz/5;
  u_NIO = @(z) U_io*exp((z/Ld));
  v_NIO = @(z) zeros(size(z));
 
  wvt.addInertialMotions(u_NIO,v_NIO);
  ```
 
  It is important to note that because the WVTransform
  de-aliases by default, you will not likely get exactly the
  same function out that you put in. The high-modes are
  removed.
 
  The new inertial motions are added to the existing inertial motions
      
