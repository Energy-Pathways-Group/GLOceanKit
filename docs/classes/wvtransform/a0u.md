---
layout: default
title: A0U
parent: WVTransform
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  A0U

matrix component that multiplies $$\tilde{u}$$ to compute $$A_0$$.


---

## Description
Complex valued transform property with dimensions $$(k,l,j)$$ and units of $$s$$.

## Discussion

These are the row 3, column 1 components of the [inverse wave-vortex (S)orting matrix](/mathematical-introduction/transformations.html), referred to as $$S^{-1}$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The primary internal gravity wave and geostrophic solutions that exist for $$k^2+l^2>0, j>0$$ are summarized in equation C5.

For $$k^2+l^2>0, j>0$$ (from either equation B14 or C5) this is written as,

$$
\textrm{A0U} \equiv i \frac{l h f_0}{\omega^2}
$$

in the manuscript. In code this is computed with,

```matlab
fOmega = f./omega;
A0U = sqrt(-1)*self.h.*(fOmega./omega) .* L;
```

With a rigid lid the solution at $$k>0, l>0, j=0$$ is from equation B11,

$$
\textrm{A0U} \equiv i \frac{f l}{g K^2}
$$

which in code is,

```matlab
A0U(:,:,1) = sqrt(-1)*(f/g_)*L(:,:,1)./K2(:,:,1);
```

The $$k=l=0, j>=0$$ solution is a mean density anomaly,

```matlab
A0U(1,1,:) = 0;
```

