---
layout: default
title: transformWaveVortexToUVWEta
parent: WVTransform
grand_parent: Classes
nav_order: 176
mathjax: true
---

#  transformWaveVortexToUVWEta

transform wave-vortex coefficients $$(A_+,A_-,A_0)$$ to fluid variables $$(u,v,\eta)$$.


---

## Declaration
```matlab
 [u,v,w,n] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)
```
## Parameters
+ `Ap`  positive wave coefficients at reference time t0
+ `Am`  negative wave coefficients at reference time t0
+ `A0`  geostrophic coefficients at reference time t0
+ `t`  (optional) time of observations

## Returns
+ `u`  x-component of the fluid velocity
+ `v`  y-component of the fluid velocity
+ `w`  z-component of the fluid velocity
+ `n`  scaled density anomaly

## Discussion

  This function is the inverse WVTransform. It is a
  [linear
  transformation](/mathematical-introduction/transformations.html)
  denoted $$\mathcal{L}$$.
 
  This function is not intended to be used directly (although
  you can), and is kept here to demonstrate a simple
  implementation of the transformation. Instead, you should
  initialize the WVTransform using one of the
  initialization functions.
 
                    
