---
layout: default
title: initWithInertialMotions
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 28
mathjax: true
---

#  initWithInertialMotions

initialize with inertial motions


---

## Declaration
```matlab
 initWithInertialMotions(self,u,v)
```
## Parameters
+ `u`  function handle that takes a single argument, u(Z)
+ `v`  function handle that takes a single argument, v(Z)

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
 
        
