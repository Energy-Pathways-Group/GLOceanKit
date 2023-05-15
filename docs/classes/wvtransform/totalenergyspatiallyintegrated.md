---
layout: default
title: totalEnergySpatiallyIntegrated
parent: WVTransform
grand_parent: Classes
nav_order: 166
mathjax: true
---

#  totalEnergySpatiallyIntegrated

horizontally-averaged depth-integrated energy computed in the spatial domain


---

## Description
Real valued state variable with no dimensions and units of $$m3/s2$$.

## Discussion
% 
The horizontally-averaged depth-integrated energy computed in the spatial domain is defined as,

$$
\frac{1}{2 L_x L_y} \int_0^{Lx} \int_0^{Ly} \int_{-L_z}^0 \left( u^2 + v^2 + w^2 + N^2 \eta^2 \right) dz
$$

In code,

```matlab
[u,v,w,eta] = self.variables('u','v','w','eta');
energy = trapz(self.z,mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
```

