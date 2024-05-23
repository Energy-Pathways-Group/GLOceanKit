---
layout: default
title: setWaveForcingCoefficients
parent: WVNonlinearFluxForced
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  setWaveForcingCoefficients

set forcing values for the wave part of the flow


---

## Declaration
```matlab
  setWaveForcingCoefficients(Apbar,Ambar,options)
```
## Parameters
+ `Apbar`  Ap 'mean' value to relax to
+ `Ambar`  Am 'mean' value to relax to
+ `MAp`  (optional) forcing mask, Ap. 1s at the forced modes, 0s at the unforced modes. Default is MAp = abs(Apbar) > 0
+ `MAm`  (optional) forcing mask, Am. 1s at the forced modes, 0s at the unforced modes. Default is MAm = abs(Apbar) > 0
+ `tauP`  (optional) relaxation time (default 0)
+ `tauM`  (optional) relaxation time (default 0)

## Discussion

  Forcing takes the following form,
 
  $$
  \frac{\partial}{\partial t} A_\pm^{klj} = \underbrace{M_{A_\pm}^{klj} \left(\bar{A}_\pm^{klj}  - A_\pm^{klj} \right)/ \tau_\pm}_{F_\textrm{force}} + F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
  $$
 
  where $$M_{A_\pm}^{klj}$$ are masks (1s and 0s),
  $$\bar{A}_\pm^{klj}$$ are mean amplitudes, and $$\tau_\pm$$
  are time scales. If the time scale is set to 0, then the mean
  amplitudes remain fixed for all time. In that limit, the
  equations can be written as,
 
  $$
  \frac{\partial}{\partial t} A_\pm^{klj} = \neg M_{A_\pm}^{klj} \left( F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj} \right)
  $$
 
                
