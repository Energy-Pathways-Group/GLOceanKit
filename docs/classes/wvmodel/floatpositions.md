---
layout: default
title: floatPositions
parent: WVModel
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  floatPositions

Returns the positions of the floats at the current time as well as the value of the fields being tracked.


---

## Declaration
```matlab
 [x,y,z,tracked] = floatPositions()
```
## Discussion

     
  The tracked variable is a structure, with fields named for
  each of the requested fields being tracked.
 
  In the following example, float positions are set along with
  one tracked field 
  ```matlab
  model.setFloatPositions(xFloat,yFloat,zFloat,'rho_total');
 
  % Set up the integrator
  nT = model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10,finalTime=3*period);
 
  % write the float trajectories to memory
  xFloatT = zeros(nT,nTrajectories);
  yFloatT = zeros(nT,nTrajectories);
  zFloatT = zeros(nT,nTrajectories);
  rhoFloatT = zeros(nT,nTrajectories);
  t = zeros(nT,1);
 
  [xFloatT(1,:),yFloatT(1,:),zFloatT(1,:),tracked] = model.floatPositions;
  rhoFloatT(1,:) = tracked.rho_total;
  ```
