---
layout: default
title: transformToKLAxes
parent: WVTransform
grand_parent: Classes
nav_order: 176
mathjax: true
---

#  transformToKLAxes

transforms in the spectral domain from (j,kl) to (kAxis,lAxis,j)


---

## Declaration
```matlab
 [varargout] = transformToKLAxes(varargin) 
```
## Parameters
+ `varargin`  variables with dimensions $$(j,kl)$$

## Returns
+ `varargout`  variables with dimensions (kAxis,lAxis,j)

## Discussion

 
  The following example takes the total energy of the wave part of
  flow, transforms it to the (kAxis,lAxis,j) grid, then sums all the energy
  along the mode (j) dimension.
 
  ```matlab
  figure
  TE = wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2);
  pcolor(wvt.kAxis,wvt.lAxis,log10(sum(wvt.transformToKLAxes(TE),3)).'), shading flat
  ```
 
        
