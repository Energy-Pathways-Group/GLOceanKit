---
layout: default
title: transformToRadialWavenumber
parent: WVTransform
grand_parent: Classes
nav_order: 172
mathjax: true
---

#  transformToRadialWavenumber

transforms in the spectral domain from (k,l,j) to (kRadial,j)


---

## Declaration
```matlab
 [varargout] = transformToRadialWavenumber(varargin) 
```
## Parameters
+ `varargin`  variables with dimensions $$(k,l)$$ or $$(k,l,j)$$

## Returns
+ `varargout`  variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

## Discussion

  Sums all the variance/energy in radial bins `kRadial`.
 
  The following example takes the total energy of the geostrophic part of
  flow, converts it to a one-dimensional spectrum in $$k$$, and then plots
  it with pcolor. The next plot then sums over all wavenumber, and produces
  plots the total energy spectrum as a function of vertical mode $$j$$
  only.
 
  ```matlab
  figure
  tiledlayout('flow')
  Ekj = wvt.transformToRadialWavenumber( wvt.A0_TE_factor .* abs(wvt.A0).^2 );
  nexttile, pcolor(wvt.j,wvt.kRadial,Ekj), shading flat
  nexttile, plot(wvt.j,sum(Ekj,1))
  ```
 
        
