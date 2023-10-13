---
layout: default
title: setFloatPositions
parent: WVModel
grand_parent: Classes
nav_order: 21
mathjax: true
---

#  setFloatPositions

Set positions of float-like particles to be advected by the model.


---

## Declaration
```matlab
 setFloatPositions(self,x,y,z,trackedFields,options)
```
## Parameters
+ `x`  x-coordinate location of the particles
+ `y`  y-coordinate location of the particles
+ `z`  z-coordinate location of the particles
+ `trackedFields`  strings of variable names
+ `advectionInterpolation`  (optional) interpolation method used for particle advection. "linear" (default), "spline", "exact"
+ `trackedVarInterpolation`  (optional) interpolation method used for tracked field. "linear" (default), "spline", "exact"

## Discussion

                 
  Pass the initial positions of particles to be advected by all
  three components of the velocity field, (u,v,w).
 
  Particles move between grid (collocation) points and thus
  their location must be interpolated. By default the
  advectionInterpolation is set to "linear" interpolation. For
  many flows this will have sufficient accuracy and allow you
  to place float at nearly every grid point without slowing
  down the model integration. However, if high accuracy is
  required, you may want to use cubic "spline" interpolation or
  even "exact" at the expense of computational speed.
 
  You can track the value of any known WVVariableAnnotation along the
  particle's flow path, e.g., relative vorticity. These values
  must also be interpolated using one of the known
  interpolation methods.
 
  ```matlab
  nTrajectories = 101;
  xFloat = Lx/2*ones(1,nTrajectories);
  yFloat = Ly/2*ones(1,nTrajectories);
  zFloat = linspace(-Lz,0,nTrajectories);
 
  model.setFloatPositions(xFloat,yFloat,zFloat,'rho_total');
  ```
 
  If a NetCDF file is set for output, the particle positions
  and tracked fields will automatically be written to file
  during integration. If you are not writing to file you can
  retrieve the current positions and values of the tracked
  fields by calling -floatPositions.
