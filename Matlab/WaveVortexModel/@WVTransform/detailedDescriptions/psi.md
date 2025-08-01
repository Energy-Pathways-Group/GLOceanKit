- topic: State Variables

The geostrophic streamfunction $$\psi$$ is computed from,

$$
\hat{\psi} = \frac{g}{f_0} A_0
$$

and then transformed back to the spatial domain with the $$F$$ modes using [transformToSpatialDomainWithF](classes/wvtransform/transformtospatialdomainwithf.html).

In code,

```matlab
f = @(wvt) wvt.transformToSpatialDomainWithF(A0=(wvt.g/wvt.f) * wvt.A0t);
```
