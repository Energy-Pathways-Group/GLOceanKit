---
layout: default
title: addInertialMotions
parent: WVTransform
grand_parent: Classes
nav_order: 61
mathjax: true
---

#  addInertialMotions

add inertial motions to existing inertial motions


---

## Declaration
```matlab
 addInertialMotions(self,u,v)
```
## Parameters
+ `u`  function handle that takes a single argument, u(Z)
+ `v`  function handle that takes a single argument, v(Z)

## Discussion

  ```matlab
  U_io = 0.2;
  Ld = wvt.Lz/5;
  u_NIO = @(z) U_io*exp(-(z/Ld));
  v_NIO = @(z) zeros(size(z));
  
  wvt.addInertialMotions(u_NIO,v_NIO);
  ```
 
  The new inertial motions are added to the existing inertial motions
        
