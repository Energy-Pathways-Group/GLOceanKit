---
layout: default
title: maskForNyquistModes
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 33
mathjax: true
---

#  maskForNyquistModes

returns a mask with locations of modes that are not fully resolved


---

## Declaration
```matlab
 nyquistMask = WVGeometryDoublyPeriodic.maskForNyquistModes(Nx,Ny,Nz);
```
## Parameters
+ `Nx`  grid points in the x-direction
+ `Ny`  grid points in the y-direction
+ `Nz`  grid points in the z-direction (defuault 1)

## Returns
+ `nyquistMask`  mask aliased mode

## Discussion

  Returns a 'mask' (matrices with 1s or 0s) indicating where Nyquist
  modes are located a standard FFT matrix.
 
  Basic usage,
  ```matlab
  NyquistMask = wvm.maskForNyquistModes(8,8);
  ```
  will return a mask that contains 1 at the locations of modes that will
  are at the Nyquist frequency of the Fourier transforms.
 
            
