- Topic: Wave-vortex sorting matrix — inverse components ($$S^{-1}$$)

These are the row 1, column 1 components of the [inverse wave-vortex (S)orting matrix](/mathematical-introduction/transformations.html), referred to as $$S^{-1}$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The primary internal gravity wave and geostrophic solutions that exist for $$k^2+l^2>0, j>0$$ are summarized in equation C5.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{ApU} \equiv \frac{k \omega + i l f_0}{2 \omega K}
$$

in the manuscript. In code this is computed with,

```matlab
alpha = atan2(L,K);
fOmega = f./omega;
ApU = (1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
ApU(:,:,1) = 0;
```

The inertial solutions occupy the $$k^2+l^2=0$$ portion of the matrix,

```matlab
ApU(1,1,:) = 1/2;
```
