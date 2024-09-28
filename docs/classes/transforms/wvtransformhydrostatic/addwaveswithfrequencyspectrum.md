---
layout: default
title: addWavesWithFrequencySpectrum
parent: WVTransformHydrostatic
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  addWavesWithFrequencySpectrum

add waves with a specified frequency spectrum


---

## Declaration
```matlab
 addWavesWithFrequencySpectrum(options)
```
## Parameters
+ `ApmSpectrum`  function_handle with signature @(omega,j), defaults to a white spectrum.
+ `shouldOnlyRandomizeOrientations`  boolean indicating whether randomness in amplitudes should be eliminated (default 0)
+ `shouldShowDiagnostics`  whether to summarize what just happened (default 0)

## Discussion

  This allows you to initialize the wave field (Ap,Am matrices)
  with a spectrum specified in terms of vertical mode j and
  frequency $\omega$. This allows us to initialize with a
  Garrett-Munk spectrum, for example, using code like,
 
  ```matlab
  GM = @(omega,j) E*H(j) .* B(omega);
  ```
  
  Because the model has limited resolution, there will not
  necessarily be many modes in a given frequency band. This
  means that the ensemble may be over a very low number of
  realization, and thus might not converge to the requested
  spectrum. For this reason, the option
  shouldOnlyRandomizeOrientations may be useful. This will only
  randomize the phases of the waves, while fixing the
  amplitudes so that the desired spectrum will be achieved.
 
          
