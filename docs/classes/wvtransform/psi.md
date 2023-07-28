---
layout: default
title: psi
parent: WVTransform
grand_parent: Classes
nav_order: 138
mathjax: true
---

#  psi

geostrophic streamfunction


---

## Description
Real valued state variable with dimensions $$(x,y,z)$$ and units of $$m^2/s$$.

## Discussion

The geostrophic streamfunction $$\psi$$ is computed from,

$$
\hat{\psi} = \frac{g}{f_0} A_0
$$

and then transformed back to the spatial domain with the $$F$$ modes using [transformToSpatialDomainWithF](classes/wvtransform/transformtospatialdomainwithf.html).

In code,

```matlab
f = @(wvt) wvt.transformToSpatialDomainWithF((wvt.g/wvt.f) * wvt.A0t);
```

