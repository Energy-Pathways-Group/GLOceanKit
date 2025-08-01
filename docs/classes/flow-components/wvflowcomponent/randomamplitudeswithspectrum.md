---
layout: default
title: randomAmplitudesWithSpectrum
parent: WVFlowComponent
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  randomAmplitudesWithSpectrum

initialize with coefficients following a specified spectrum


---

## Declaration
```matlab
 Ap,Am,A0] = randomAmplitudesWithSpectrum(options)
```
## Parameters
+ `A0Spectrum`  (optional) function_handle with signature @(k,j), defaults to a white spectrum.
+ `ApmSpectrum`  (optional) function_handle with signature @(k,j), defaults to a white spectrum.
+ `shouldOnlyRandomizeOrientations`  boolean indicating whether randomness in amplitudes should be eliminated (default 0)

## Returns
+ `Ap`  matrix of size [Nj Nkl]
+ `Am`  matrix of size [Nj Nkl]
+ `A0`  matrix of size [Nj Nkl]

## Discussion

  This allows you to initialize amplitudes following a spectrum
  defined in terms of wavenumber and vertical mode.
 
                
