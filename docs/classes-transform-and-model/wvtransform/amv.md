---
layout: default
title: AmV
parent: WVTransform
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  AmV

matrix component that multiplies $$\tilde{v}$$ to compute $$A_m$$.


---

## Description
Complex valued transform property with dimensions $$(k,l,j)$$ and no units.

## Discussion

These are the row 2, column 2 components of the [inverse wave-vortex (S)orting matrix](/mathematical-introduction/transformations.html), referred to as $$S^{-1}$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The primary internal gravity wave and geostrophic solutions that exist for $$k^2+l^2>0, j>0$$ are summarized in equation C5.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{AmV} \equiv \frac{l \omega + i k f_0}{2 \omega K}
$$

in the manuscript. In code this is computed with,

```matlab
alpha = atan2(L,K);
fOmega = f./omega;
AmV = (1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
AmV(:,:,1) = 0;
```

The inertial solutions occupy the $$k^2+l^2=0$$ portion of the matrix,

```matlab
AmV(1,1,:) = sqrt(-1)/2;
```

