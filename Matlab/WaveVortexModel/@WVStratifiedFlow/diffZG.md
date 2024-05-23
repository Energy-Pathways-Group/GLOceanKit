- Topic: Operations — Differentiation
- Declaration: wz = diffZG(w)
- Parameter w: variable with dimensions $$(x,y,z)$$
- Returns wz: differentiated variable with dimensions $$(x,y,z)$$

Each subclass implements this operation differently, depending on the vertical modes being used.

For hydrostatic vertical modes with a rigid-lid and zero buoyancy anomaly,

$$
\partial_z w = \mathcal{F}^{-1} \left[ \frac{1}{h_j} \mathcal{G}\left[ w \right] \right]
$$

where we've used the same notation as defined for the [discrete transformations](/mathematical-introduction/transformations.html).
