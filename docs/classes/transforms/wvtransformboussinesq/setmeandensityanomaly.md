---
layout: default
title: setMeanDensityAnomaly
parent: WVTransformBoussinesq
grand_parent: Classes
nav_order: 44
mathjax: true
---

#  setMeanDensityAnomaly

set inertial motions


---

## Declaration
```matlab
 setMeanDensityAnomaly(eta)
```
## Parameters
+ `eta`  function handle that takes a single argument, eta(Z)

## Discussion

  Overwrites existing inertial motions with the new values.
  Other components of the flow will remain unaffected.
 
  ```matlab
  U_io = 0.2;
  Ld = wvt.Lz/5;
  u_NIO = @(z) U_io*exp((z/Ld));
  v_NIO = @(z) zeros(size(z));
 
  wvt.setInertialMotions(u_NIO,v_NIO);
  ```
 
  It is important to note that because the WVTransform
  de-aliases by default, you will not likely get exactly the
  same function out that you put in. The high-modes are
  removed.
             
      
