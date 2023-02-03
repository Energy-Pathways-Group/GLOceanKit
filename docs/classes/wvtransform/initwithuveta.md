---
layout: default
title: initWithUVEta
parent: WVTransform
grand_parent: Classes
nav_order: 107
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

  Replaces the variables Ap,Am,A0 with those computed from $$(u,v,\eta)$$.
  If a time t is specified, the wvt is set to that time.
            
