---
layout: default
title: spatialFlux
parent: WVNonlinearFlux
grand_parent: Classes
nav_order: 17
mathjax: true
---

#  spatialFlux

a subclass can override this, and then modify the spatial


---

## Discussion
fluxes that get returned.
 
  To count everything:
  1. (u,v,w,eta) -> 4 transformToSpatialDomainWithFourier
  2. (u_x,u_y,v_x,v_y,eta_x,eta_y) -> 3 diffX, 3 diffY
  3. (Ap,Am,A0) -> 3 transformFomSpatialDomainWithFourier
