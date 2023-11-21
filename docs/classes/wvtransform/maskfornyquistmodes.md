---
layout: default
title: maskForNyquistModes
parent: WVTransform
grand_parent: Classes
nav_order: 132
mathjax: true
---

#  maskForNyquistModes

returns a mask with locations of modes that are not fully resolved


---

## Declaration
```matlab
 NyquistMask = maskForNyquistModes();
```
## Returns
+ `NyquistMask`  mask nyquist modes

## Discussion

  Returns a 'mask' (matrices with 1s or 0s) indicating where Nyquist 
  modes are located in the Ap/Am/A0 matrices.
 
  Basic usage,
  NyquistMask = wvm.maskForNyquistModes();
  will return a mask that contains 1 at the locations of modes that will
  are at the Nyquist frequency of the Fourier transforms.
 
      
