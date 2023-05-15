---
layout: default
title: transformToRadialWavenumber
parent: WVTransform
grand_parent: Classes
nav_order: 170
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
 
        
