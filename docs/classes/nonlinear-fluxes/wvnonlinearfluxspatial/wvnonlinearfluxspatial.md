---
layout: default
title: WVNonlinearFluxSpatial
parent: WVNonlinearFluxSpatial
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  WVNonlinearFluxSpatial

3D nonlinear flux for Boussinesq flow, computed in the spatial domain


---

## Declaration
```matlab
 WVNonlinearFluxSpatial < [WVNonlinearFluxOperation](/classes/wvnonlinearfluxoperation/)
```
## Discussion

  Computes the nonlinear flux for a Boussinesq model. This class is not
  intended to be used for numerical modeling as it does not have any
  antialiasing or damping, but is indended as an example. The
  implementation is *simple* and follows directly from the equations of
  motion, but it is not the fastest implementation. To compute
  nonlinear fluxes appropriate for numerical modeling, use the
  [WVNonlinearFlux](/classes/wvnonlinearflux/) class.
 
    
