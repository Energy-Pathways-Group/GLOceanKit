- Topic: Operations — Differentiation
- Declaration: uz = diffZF(u)
- Parameter u: variable with dimensions $$(x,y,z)$$
- Returns uz: differentiated variable with dimensions $$(x,y,z)$$

Each subclass implements this operation differently, depending on the vertical modes being used.

For hydrostatic vertical modes with a rigid-lid and zero buoyancy anomaly,

$$
\partial_z u = -\frac{1}{g} N^2(z)  \mathcal{G}^{-1} \left[ \mathcal{F} \left[ u \right] \right]
$$

where we've used the same notation as defined for the [discrete transformations](/mathematical-introduction/transformations.html).
