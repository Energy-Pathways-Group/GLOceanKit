---
layout: default
title: setInertialMotions
parent: WVTransform
grand_parent: Classes
nav_order: 157
mathjax: true
---

#  setInertialMotions

set inertial motions


---

## Declaration
```matlab
 setInertialMotions(self,u,v)
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
  
  wvt.setInertialMotions(u_NIO,v_NIO);
  ```
 
  Overwrites existing inertial motions with the new values
        
