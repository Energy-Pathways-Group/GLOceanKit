- Topic: Operations — Transformations
- Declaration: [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(u_bar)
- Parameter u_bar: variable with dimensions $$(k,l,j)$$
- Returns u: variable u with dimensions $$(x,y,z)$$
- Returns ux: variable du/dx with dimensions $$(x,y,z)$$
- Returns uy: variable du/dy with dimensions $$(x,y,z)$$
- Returns uz: variable du/dz with dimensions $$(x,y,z)$$

This performs the same operation as `transformToSpatialDomainWithF`, but also returns the first-derivative in all three spatial directions.

The computation of these derivatives can be performed more efficiently if done simultaneously. So when performance is a requirement, it can be useful to call this function rather than request the derivatives individually.
