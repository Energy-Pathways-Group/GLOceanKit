- Topic: Wave-vortex sorting matrix — inverse components ($$S^{-1}$$)

These are the row 3, column 3 components of the [inverse wave-vortex (S)orting matrix](/mathematical-introduction/transformations.html), referred to as $$S^{-1}$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The primary internal gravity wave and geostrophic solutions that exist for $$k^2+l^2>0, j>0$$ are summarized in equation C5.

For $$k^2+l^2>0, j>0$$ (from either equation B14 or C5) this is written as,

$$
\textrm{A0N} \equiv \frac{f_0^2}{\omega^2}
$$

in the manuscript. In code this is computed with,

```matlab
fOmega = f./omega;
A0N = fOmega.^2;
```

With a rigid lid the solution at $$k>0, l>0, j=0$$ is from equation B11,

$$
\textrm{A0N} \equiv 0
$$

which in code is,

```matlab
A0N(:,:,1) = 0;
```

The $$k=l=0, j>=0$$ solution is a mean density anomaly,

```matlab
A0N(1,1,:) = 1;
```
