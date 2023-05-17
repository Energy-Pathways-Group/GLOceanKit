---
layout: default
title: initWithSpectralFunction
parent: WVTransform
grand_parent: Classes
nav_order: 106
mathjax: true
---

#  initWithSpectralFunction

initialize the wave spectrum with a given function


---

## Declaration
```matlab
 [GM3Dint,GM3Dext] = initWithSpectralFunction(GM2D_int, varargin) 
```
## Discussion

     
  The GM2D_int function is used to assign variance to a given
  wave mode. It has three arguments, omega0, omega1, and j and
  should return the amount of variance you want assigned to a
  wave mode between omega0 and omega1 at vertical mode j.
 
  The returned values GM3Dint are the results of distributing
  this variance. size(GM3Dint) = size(Kh), so you can see
  how much energy was assigned to each internal mode and
  similarly size(GM3Dext) = size(self.k_ext).
 
  The function takes the (optional) name/value pairs:
 
  shouldRandomizeAmplitude = 1 or 0 will randomize the
  energy in each mode such that the expected value matches that
  assigned. Default 0 (amplitudes will not be randomized)
 
  maxDeltaOmega is the maximum width in frequency that will be
  integrated over for assigned energy. By default it is self.Nmax-self.f
