---
layout: default
title: initWithUVRho
parent: WVTransform
grand_parent: Classes
nav_order: 110
mathjax: true
---

#  initWithUVRho

initialize with fluid variables $$(u,v,\rho)$$


---

## Declaration
```matlab
 initWithUVRho(U,V,N,t)
```
## Parameters
+ `u`  x-component of the fluid velocity
+ `v`  y-component of the fluid velocity
+ `n`  scaled density anomaly
+ `t`  (optional) time of observations

## Discussion

  Clears variables Ap,Am,A0 and then randomizes the flow
            
