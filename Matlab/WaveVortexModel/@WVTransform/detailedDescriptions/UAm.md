- Topic: Wave-vortex sorting matrix — components of $$S$$

These are the row 1, column 2 components of the [wave-vortex (S)orting matrix](/transformations/transformations.html), referred to as the $$S$$ matrix in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The primary internal gravity wave and geostrophic solutions that exist for $$k^2+l^2>0, j>0$$ are summarized in equation C4.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{UAm} \equiv \frac{k \omega + i l f_0}{\omega K}
$$

in the manuscript. In code this is computed with,

```matlab
alpha = atan2(L,K);
fOmega = f./omega;
UAm = (cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
UAm(:,:,1) = 0;
```

The inertial solutions occupy the $$k^2+l^2=0$$ portion of the matrix,

```matlab
UAm(1,1,:) = 1;
```