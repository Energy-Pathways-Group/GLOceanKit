---
layout: default
title: setupIntegrator
parent: WVModel
grand_parent: Classes
nav_order: 40
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
+ `integratorType`  (optional) integrator type. "adaptive"(default), "fixed"
+ `deltaT`  (fixed) time step
+ `cfl`  (fixed) cfl condition
+ `timeStepConstraint`  (fixed) constraint to fix the time step. "advective" (default) ,"oscillatory","min"
+ `integrator`  (adapative) function handle of integrator. @ode78 (default)
+ `absTolerance`  (adapative) absolute tolerance for sqrt(energy). 1e-6 (default)
+ `relTolerance`  (adapative) relative tolerance for sqrt(energy). 1e-3 (default)
+ `shouldUseScaledTolerance`  (adapative) whether to scale by the energy norm. 0 or 1 (default)
+ `absToleranceA0`  (adapative) absolute tolerance for A0 used when shouldUseScaledTolerance=0 . 1e-10 (default)
+ `absToleranceApm`  (adapative) absolute tolerance for Apm used when shouldUseScaledTolerance=0 . 1e-10 (default)
+ `absToleranceXY`  (adapative) absolute tolerance in meters for particle advection in (x,y). 1e-1 (default)
+ `absToleranceZ`  (adapative) absolute tolerance  in meters for particle advection in (z). 1e-2 (default)
+ `shouldShowIntegrationStats`  (adapative) whether to show integration output 0 or 1 (default)

## Discussion

  By default the model will use adaptive time stepping with a
  reasonable choice of values. However, you may find it
  necessary to customize the time stepping behavior.
 
  When setting up the integrator you must choice between
  "adaptive" and "fixed" integrator types. Depending on which
  type you choose, you will have different options available.
 
  The "fixed" time-step integrator used a cfl condition based
  on the advective velocity, but you can change this to use the
  highest oscillatory frequency. Alternatively, you can simply
  set deltaT yourself.
 
  The "adaptive" time-step integator uses absolute and relative
  error tolerances. It is worth reading Matlab's documentation
  on RelTol and AbsTol as part of odeset to understand what
  these mean. By default, the adaptive time stepping uses a
  a relative error tolerance of 1e-3 for everything. However,
  the absolute error tolerance is less straightforward.
 
  The absolute tolerance has a meaningful scale with units, and
  thus must be chosen differently for particle positions (x,y)
  than for geostrophic coefficients (A0). 
 
                              
