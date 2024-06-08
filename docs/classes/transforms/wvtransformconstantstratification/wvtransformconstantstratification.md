---
layout: default
title: WVTransformConstantStratification
parent: WVTransformConstantStratification
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  WVTransformConstantStratification

initialze a wave-vortex transform with constant stratification


---

## Declaration
```matlab
 wvt = WVTransformConstantStratification(Lxyz, Nxyz, options)
```
## Parameters
+ `Lxyz`  length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
+ `Nxyz`  number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
+ `N0`   (optional) buoyancy frequency (radians/s) default is 5.2e-3, or 3 cph)
+ `latitude`  (optional) latitude of the domain (default is 33 degrees north)
+ `rho0`  (optional) density at the surface z=0 (default is 1025 kg/m^3)
+ `isHydrostatic`  (optional) flag indicating whether to use hydrostatic transformations (default 0)

## Returns
+ `wvt`  a new WVTransformConstantStratification instance

## Discussion

 
  To initialization an instance of the
  WVTransformConstantStratification class you must specific the
  stratification and latitude, otherwise defaults will be
  assumed for you.
 
  To initialize a domain of size (100 km, 100 km, 1300 m) with
  a buoyancy frequency of 3 cycles per hour at latitude 30,
  call
 
  ```matlab
  N0 = 3*2*pi/3600;
  wvt = WVTransformConstantStratification([100e3, 100e3, 1300],[64, 64, 65], N0=N0,latitude=30);
  ```
 
 
 
                  
