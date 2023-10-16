---
layout: default
title: UAp
parent: WVTransform
grand_parent: Classes
nav_order: 48
mathjax: true
---

#  UAp

matrix component that multiplies $$A_p$$ to compute $$\tilde{u}$$.


---

## Description
Complex valued transform property with dimensions $$(k,l,j)$$ and no units.

## Discussion

These are the row 1, column 1 components of the [wave-vortex (S)orting matrix](/mathematical-introduction/transformations.html), referred to as the $$S$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The primary internal gravity wave and geostrophic solutions that exist for $$k^2+l^2>0, j>0$$ are summarized in equation C4.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{UAp} \equiv \frac{k \omega - i l f_0}{\omega K}
$$

in the manuscript. In code this is computed with,

```matlab
alpha = atan2(L,K);
fOmega = f./omega;
UAp = (cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
UAp(:,:,1) = 0;
```

The inertial solutions occupy the $$k^2+l^2=0$$ portion of the matrix,

```matlab
UAp(1,1,:) = 1;
```

