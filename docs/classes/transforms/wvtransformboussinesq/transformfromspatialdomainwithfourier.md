---
layout: default
title: transformFromSpatialDomainWithFourier
parent: WVTransformBoussinesq
grand_parent: Classes
nav_order: 50
mathjax: true
---

#  transformFromSpatialDomainWithFourier

self.dftBuffer = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);


---

## Discussion
u_bar(self.wvPrimaryIndex) = self.dftBuffer(self.dftPrimaryIndex);
  u_bar=reshape(u_bar,[self.Nz self.Nkl]);
