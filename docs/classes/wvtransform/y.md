---
layout: default
title: y
parent: WVTransform
grand_parent: Classes
nav_order: 224
mathjax: true
---

#  y

y coordinate


---

## Discussion

The values `Ly` and `Ny` are set during initialization from which the `y` coordinate is derived.

The y coordinate is periodic, which means that
```matlab
dy = Ly/Ny;
y = dy*(0:Ny-1)';
```

Note that this means that it is NOT true that Ly=y(end)-y(1), but in fact you need an extra grid point, i.e., it IS true that Ly = dy + y(end)-y(1). This is the usual grid for Fourier Transforms.

