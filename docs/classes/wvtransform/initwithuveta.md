---
layout: default
title: initWithUVEta
parent: WVTransform
grand_parent: Classes
nav_order: 108
mathjax: true
---

#  initWithUVEta

initialize with fluid variables $$(u,v,\eta)$$


---

## Declaration
```matlab
 initWithUVEta(U,V,N,t)
```
## Parameters
+ `u`  x-component of the fluid velocity
+ `v`  y-component of the fluid velocity
+ `n`  scaled density anomaly
+ `t`  (optional) time of observations

## Discussion

  Clears variables Ap,Am,A0 and then randomizes the flow
            