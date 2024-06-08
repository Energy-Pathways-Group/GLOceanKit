---
layout: default
title: initWithMeanDensityAnomaly
parent: WVTransformBoussinesq
grand_parent: Classes
nav_order: 28
mathjax: true
---

#  initWithMeanDensityAnomaly

initialize with inertial motions


---

## Declaration
```matlab
 initWithMeanDensityAnomaly(eta)
```
## Parameters
+ `eta`  function handle that takes a single argument, eta(Z)

## Discussion

  Clears variables Ap,Am,A0 and then sets inertial motions.
  
  ```matlab
  U_io = 0.2;
  Ld = wvt.Lz/5;
  u_NIO = @(z) U_io*exp((z/Ld));
  v_NIO = @(z) zeros(size(z));
 
  wvt.initWithInertialMotions(u_NIO,v_NIO);
  ```
 
  It is important to note that because the WVTransform
  de-aliases by default, you will not likely get exactly the
  same function out that you put in. The high-modes are
  removed.
 
      
