---
layout: default
title: setupIntegrator
parent: WVModel
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  setupIntegrator

Customize the time-stepping


---

## Declaration
```matlab
 setupIntegrator(self,options)
```
## Parameters
+ `deltaT`  (optional) set the integrator time step
+ `cfl`  (optional) set the cfl condition used to set integrator time step
+ `timeStepConstraint`  (optional) set the method used to determine the integrator time step. "advective","oscillatory","min"
+ `outputInterval`  (optional) If set, it will allow you to call -integrateToNextOutputTime and, if a NetCDF file is set for output, it will set the interval at which time steps are written to file.
+ `finalTime`  (optional) if set, the NetCDF file may set a fixed time dimension length.

## Discussion

              
