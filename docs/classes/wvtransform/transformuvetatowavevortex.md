---
layout: default
title: transformUVEtaToWaveVortex
parent: WVTransform
grand_parent: Classes
nav_order: 177
mathjax: true
---

#  transformUVEtaToWaveVortex

transform fluid variables $$(u,v,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.


---

## Declaration
```matlab
 [Ap,Am,A0] = transformUVEtaToWaveVortex(U,V,N,t)
```
## Parameters
+ `u`  x-component of the fluid velocity
+ `v`  y-component of the fluid velocity
+ `n`  scaled density anomaly
+ `t`  (optional) time of observations

## Returns
+ `Ap`  positive wave coefficients at reference time t0
+ `Am`  negative wave coefficients at reference time t0
+ `A0`  geostrophic coefficients at reference time t0

## Discussion

  This function **is** the WVTransform. It is a [linear
  transformation](/mathematical-introduction/transformations.html)
  denoted $$\mathcal{L}$$.
 
  This function is not intended to be used directly (although
  you can), and is kept here to demonstrate a simple
  implementation of the transformation. Instead, you should
  initialize the WVTransform using one of the
  initialization functions.
 
                  
