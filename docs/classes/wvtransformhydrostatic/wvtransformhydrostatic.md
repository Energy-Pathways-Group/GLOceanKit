---
layout: default
title: WVTransformHydrostatic
parent: WVTransformHydrostatic
grand_parent: Classes
nav_order: 19
mathjax: true
---

#  WVTransformHydrostatic

create a wave-vortex transform for variable stratification


---

## Declaration
```matlab
 wvt = WVTransformHydrostatic(Lxyz, Nxyz, options)
```
## Parameters
+ `Lxyz`  length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
+ `Nxyz`  number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
+ `rho`   (optional) function_handle specifying the density as a function of depth on the domain [-Lz 0]
+ `stratification`   (optional) function_handle specifying the stratification as a function of depth on the domain [-Lz 0]
+ `latitude`  (optional) latitude of the domain (default is 33 degrees north)
+ `rho0`  (optional) density at the surface z=0 (default is 1025 kg/m^3)

## Returns
+ `wvt`  a new WVTransformHydrostatic instance

## Discussion

  Creates a new instance of the WVTransformHydrostatic class
  appropriate for disentangling hydrostatic waves and vortices
  in variable stratification
 
  You must initialization by passing *either* the density
  profile or the stratification profile.
 
                  
