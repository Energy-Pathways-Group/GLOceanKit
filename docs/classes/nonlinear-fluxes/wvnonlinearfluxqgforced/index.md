---
layout: default
title: WVNonlinearFluxQGForced
has_children: false
has_toc: false
mathjax: true
parent: Nonlinear fluxes
grand_parent: Class documentation
nav_order: 5
---

#  WVNonlinearFluxQGForced

3D forced quasigeostrophic potential vorticity flux


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>WVNonlinearFluxQGForced < <a href="/classes/wvnonlinearfluxqg/" title="WVNonlinearFluxQG">WVNonlinearFluxQG</a></code></pre></div></div>

## Overview
 
  The 3D quasigeostrophic potential vorticity flux will only use and
  modify the A0 coefficients.
 
  $$
  \frac{\partial}{\partial t} A_0^{klj} = \underbrace{M_{A_0}^{klj} \left(\bar{A}_0^{klj}  - A_0^{klj} \right)/ \tau_0}_{F_\textrm{force}} + F_0^{klj} + F_\textrm{damp}^{klj}
  $$
 
  To initialize the WVNonlinearFluxQGForced,
 
  ```matlab
  model = WVModel(wvt,nonlinearFlux=WVNonlinearFluxQGForced(wvt,shouldUseBeta=1,uv_damp=wvt.uMax));
  ```
 
    


## Topics
+ Initializing
+ Other
  + [`A0bar`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/a0bar.html) A0 'mean' value to relax to
  + [`MA0`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/ma0.html) Forcing mask, A0. 1s at the forced modes, 0s at the unforced modes
  + [`WVNonlinearFluxQGForced`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/wvnonlinearfluxqgforced.html) initialize 3D quasigeostrophic potential vorticity flux
  + [`compute`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/compute.html) Apply operator S---defined in (C4) in the manuscript
  + [`nonlinearFluxWithResolutionForTransform`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/nonlinearfluxwithresolutionfortransform.html) 
  + [`setNarrowBandForcing`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/setnarrowbandforcing.html) 
  + [`tau0`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/tau0.html) relaxation time
+ Equality
  + [`isequal`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxFromFile`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Computation
  + [`setGeostrophicForcingCoefficients`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/setgeostrophicforcingcoefficients.html) set forcing values for the geostrophic part of the flow
+ Write to file
  + [`writeToFile`](/classes/nonlinear-fluxes/wvnonlinearfluxqgforced/writetofile.html) write information about the nonlinear flux operation to file


---