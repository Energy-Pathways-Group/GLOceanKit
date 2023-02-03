---
layout: default
title: k
parent: WVTransform
grand_parent: Classes
nav_order: 114
mathjax: true
---

#  k

wavenumber-coordinate dimension in the x-direction


---

## Discussion

This is the usual definition of wavenumbers (or frequencies) for the Fourier transform,
```matlab
dk = 1/self.Lx; 
k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
```

The negative wavenumbers follow the positive wavenumbers. Matlab's built in function `fftshift` is useful for making these monotonically increasing.

