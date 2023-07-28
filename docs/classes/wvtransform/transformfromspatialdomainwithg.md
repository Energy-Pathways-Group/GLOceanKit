---
layout: default
title: transformFromSpatialDomainWithG
parent: WVTransform
grand_parent: Classes
nav_order: 171
mathjax: true
---

#  transformFromSpatialDomainWithG

transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the G-modes


---

## Declaration
```matlab
 w_bar = transformFromSpatialDomainWithG(w)
```
## Parameters
+ `w`  variable with dimensions $$(x,y,z)$$

## Discussion

This is the component of the [discrete transformation](/mathematical-introduction/transformations.html) $$D$$ that transforms $$(x,y) \mapsto (k,l)$$ with a discrete Fourier transform, followed by a projection onto the G-modes. Mathematically we write,

$$
\tilde{f}_{klj} =  \mathcal{G} \cdot \mathcal{DFT}_y \cdot \mathcal{DFT}_x \left[ f(x,y,z) \right]
$$

The $$G$$ mode projection is applicable to dynamical variables $$\eta$$, $$w$$.

As noted in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995), the vertical transforms $\mathcal{F}$ and $\mathcal{G}$ require a matrix multiplication and thus have a computational cost of,

$$
N_z^2 N_x N_y
$$  

while the horizontal transforms use an FFT algorithm and therefore scale as,

$$
\frac{5}{2} N_z N_x N_y \log_2 N_x N_y - N_z N_x N_y
$$

Assuming that $$N_x = N_y$$, the total computational cost of the horizontal and vertical transforms are approximately equal when $$10 log_2 N_x = N_z$$ , or $$13 log_2 N_x = N_z$$ for the hydrostatic case. This means that for a horizontal resolution of $$256^2$$ the horizontal transformations dominate the computation time until approximately $$80-100$$ vertical modes are used.

