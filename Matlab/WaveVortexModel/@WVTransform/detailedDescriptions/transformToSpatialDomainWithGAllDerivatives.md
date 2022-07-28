- Topic: Operations — Transformations
- Declaration: [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives( w_bar )
- Parameter w_bar: variable with dimensions $$(k,l,j)$$
- Returns w: variable w with dimensions $$(x,y,z)$$
- Returns wx: variable dw/dx with dimensions $$(x,y,z)$$
- Returns wy: variable dw/dy with dimensions $$(x,y,z)$$
- Returns wz: variable dw/dz with dimensions $$(x,y,z)$$

This performs the same operation as `transformToSpatialDomainWithG`, but also returns the first-derivative in all three spatial directions.

The computation of these derivatives can be performed more efficiently if done simultaneously. So when performance is a requirement, it can be useful to call this function rather than request the derivatives individually.
