---
layout: default
title: WVNonlinearFluxForced
has_children: false
has_toc: false
mathjax: true
parent: Nonlinear fluxes
grand_parent: Class documentation
nav_order: 3
---

#  WVNonlinearFluxForced

3D forced nonlinear flux for Boussinesq flow


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>WVNonlinearFluxForced < <a href="/classes/wvnonlinearflux/" title="WVNonlinearFlux">WVNonlinearFlux</a></code></pre></div></div>

## Overview
 
  The unforced model basically looks likes like this,
 
  $$
  \frac{\partial}{\partial t} A^{klj} = F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
  $$
 
  for each of the three components. The forcing adds a new term,
 
  $$
  \frac{\partial}{\partial t} A^{klj} = \underbrace{M_{A}^{klj} \left(\bar{A}^{klj}  - A^{klj} \right)/ \tau}_{F_\textrm{force}} + F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
  $$
 
  which forces those select modes to relax to their $$\bar{A}^{klj}$$
  state with time scale $$\tau$$.  If the time scale is set to 0, then the mean
  amplitudes remain fixed for all time. In that limit, the
  equations can be written as,
 
  $$
  \frac{\partial}{\partial t} A^{klj} = \neg M_{A}^{klj} \left( F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj} \right)
  $$
 
  This is most often used when initializing a model, e.g.,
 
  ```matlab
  model = WVModel(wvt,nonlinearFlux=WVNonlinearFluxForced(wvt,uv_damp=wvt.uMax));
  ```
 
    


## Topics
+ Initializing
+ Other
  + [`A0bar`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/a0bar.html) A0 'mean' value to relax to
  + [`Ambar`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/ambar.html) Am 'mean' value to relax to
  + [`Apbar`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/apbar.html) Ap 'mean' value to relax to
  + [`MA0`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/ma0.html) Forcing mask, A0. 1s at the forced modes, 0s at the unforced modes
  + [`MAm`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/mam.html) Forcing mask, Am. 1s at the forced modes, 0s at the unforced modes
  + [`MAp`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/map.html) Forcing mask, Ap. 1s at the forced modes, 0s at the unforced modes
  + [`addVariableOfType`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/addvariableoftype.html) 
  + [`tau0`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/tau0.html) A0 relaxation time
  + [`tauM`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/taum.html) Am relaxation time
  + [`tauP`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/taup.html) Ap relaxation time
+ Initialization
  + [`WVNonlinearFluxForced`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/wvnonlinearfluxforced.html) initialize WVNonlinearFluxForced
  + [`nonlinearFluxFromFile`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Computation
  + [`compute`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/compute.html) the promised variable
+ Equality
  + [`isequal`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/isequal.html) check for equality with another nonlinear flux operation
+ Set forcing
  + [`setGeostrophicForcingCoefficients`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/setgeostrophicforcingcoefficients.html) set forcing values for the geostrophic part of the flow
  + [`setWaveForcingCoefficients`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/setwaveforcingcoefficients.html) set forcing values for the wave part of the flow
+ Write to file
  + [`writeToFile`](/classes/nonlinear-fluxes/wvnonlinearfluxforced/writetofile.html) write information about the nonlinear flux operation to file


---