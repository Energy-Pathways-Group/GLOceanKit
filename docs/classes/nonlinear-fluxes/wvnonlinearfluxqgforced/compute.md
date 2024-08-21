---
layout: default
title: compute
parent: WVNonlinearFluxQGForced
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  compute

this is ever so slightly faster (for barotropic only), but why add the complication?


---

## Discussion
PVbar = self.PVA0 .* wvt.A0;
  PVx = wvt.transformToSpatialDomainWithF(A0=sqrt(-1)*shiftdim(wvt.k,-1).*PVbar);
  PVy = wvt.transformToSpatialDomainWithF(A0=sqrt(-1)*shiftdim(wvt.l,-1).*PVbar);

Help for WVNonlinearFluxQGForced/compute is inherited from superclass WVNonlinearFluxQG
