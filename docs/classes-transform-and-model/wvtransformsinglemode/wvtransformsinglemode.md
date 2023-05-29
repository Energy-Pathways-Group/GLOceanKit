---
layout: default
title: WVTransformSingleMode
parent: WVTransformSingleMode
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  WVTransformSingleMode

create a single mode wave-vortex transform


---

## Declaration
```matlab
 wvt = WVTransformSingleMode(Lxyz, Nxyz, options)
```
## Parameters
+ `Lxy`  length of the domain (in meters) in the two coordinate directions, e.g. [Lx Ly]
+ `Nxy`  number of grid points in the two coordinate directions, e.g. [Nx Ny]
+ `h`   (optional) equivalent depth (default 0.8)
+ `latitude`  (optional) latitude of the domain (default is 33 degrees north)

## Returns
+ `wvt`  a new WVTransformSingleMode instance

## Discussion

  ```matlab
  Lxy = 50e3;
  Nxy = 256;
  latitude = 25;
  wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);
  ```
 
 
              
