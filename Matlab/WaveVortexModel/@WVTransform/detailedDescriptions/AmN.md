- Topic: Wave-vortex sorting matrix — inverse components ($$S^{-1}$$)

These are the row 1, column 3 components of the [inverse wave-vortex (S)orting matrix](/mathematical-introduction/transformations.html), referred to as $$S^{-1}$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The primary internal gravity wave and geostrophic solutions that exist for $$k^2+l^2>0, j>0$$ are summarized in equation C5.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{ApN} \equiv - \frac{g K}{2 \omega}
$$

in the manuscript. In code this is computed with,

```matlab
Kh = sqrt(K.*K + L.*L);
ApN = -g*Kh./(2*omega);
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
ApN(:,:,1) = 0;
```

The inertial solutions at $$k^2+l^2=0$$ do not contribute to $$\eta$$, so that component remains zero.
