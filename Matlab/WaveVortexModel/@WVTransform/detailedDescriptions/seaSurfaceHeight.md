- topic: State Variables

The sea-surface height $$\zeta$$ for the rigid-lid model is computed by assuming that it is proportional to the pressure at the surface, $$p(z=0) = \rho_0 g \zeta$$ 

$$
\hat{\zeta} = - \frac{K h}{\omega} A_p + \frac{K h}{\omega} A_m + A_0 
$$

which is then transformed back in the spatial domain using the $$F$$ inverse transformation, [transformToSpatialDomainWithF](classes/wvtransform/transformtospatialdomainwithf.html).
