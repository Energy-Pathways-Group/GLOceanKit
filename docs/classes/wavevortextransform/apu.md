---
layout: default
title: ApU
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 18
mathjax: true
---

#  ApU

matrix coefficient that multiplies $$\tilde{u}$$ to compute $$A_p$$.


---

## Description
Real valued transform property with dimensions $$(k,l,j)$$ and no units.

## Discussion

These are the coefficients of row 1, column 1 of the $$S^{-1}$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995).

In the manuscript this is written as,

$$
\textrm{ApU} \equiv \frac{k \omega + i l f_0}{2 \omega K}
$$

and in code, this is computed with,

```matlab
alpha = atan2(L,K);
fOmega = f./omega;
ApU = (1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
```

but there are no $$k^2+l^2>0, j=0$$ wave solutions

```matlab
ApU(:,:,1) = 0;
```

inertial solutions,

```matlab
ApU(1,1,:) = 1/2;
```

