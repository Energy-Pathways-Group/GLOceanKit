---
layout: default
title: addParticles
parent: WaveVortexModel
grand_parent: Classes
nav_order: 17
mathjax: true
---

#  addParticles

Add particles to be advected by the flow.


---

## Declaration
```matlab
 addParticles(name,fluxOp,x,y,z,trackedFieldNames,options)
```
## Parameters
+ `name`  a unique name to call the particles
+ `fluxOp`  a ParticleFluxOperation, used to determine how the flow advects the particles
+ `x`  x-coordinate location of the particles
+ `y`  y-coordinate location of the particles
+ `z`  z-coordinate location of the particles
+ `trackedFields`  strings of variable names
+ `advectionInterpolation`  (optional) interpolation method used for particle advection. "linear" (default), "spline", "exact"
+ `trackedVarInterpolation`  (optional) interpolation method used for tracked field. "linear" (default), "spline", "exact"

## Discussion

                    
