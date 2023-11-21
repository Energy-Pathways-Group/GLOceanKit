---
layout: default
title: WVNonlinearFluxQG
has_children: false
has_toc: false
mathjax: true
parent: Nonlinear fluxes
grand_parent: Class documentation
nav_order: 4
---

#  WVNonlinearFluxQG

3D quasigeostrophic potential vorticity flux


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>WVNonlinearFluxQG < <a href="/classes/wvnonlinearfluxoperation/" title="WVNonlinearFluxOperation">WVNonlinearFluxOperation</a></code></pre></div></div>

## Overview
 
  The 3D quasigeostrophic potential vorticity flux will only use and
  modify the A0 coefficients.
 
  To initialize the WVNonlinearFluxQG,
 
  ```matlab
  model = WVModel(wvt,nonlinearFlux=WVNonlinearFluxQG(wvt,shouldUseBeta=1,uv_damp=wvt.uMax));
  ```
 
    


## Topics
+ Initializing
+ Other
  + [`A0PV`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/a0pv.html) conversion from PV to A0
  + [`PVA0`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/pva0.html) conversion from A0 to PV
  + [`RVA0`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/rva0.html) conversion from A0 to RV
  + [`WVNonlinearFluxQG`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/wvnonlinearfluxqg.html) initialize 3D quasigeostrophic potential vorticity flux
  + [`beta`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/beta.html) 
  + [`buildDampingOperator`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/builddampingoperator.html) 
  + [`compute`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/compute.html) Apply operator S---defined in (C4) in the manuscript
  + [`damp`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/damp.html) 
  + [`dampingOperator`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/dampingoperator.html) 
  + [`dampingTimeScale`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/dampingtimescale.html) 
  + [`nu_xy`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/nu_xy.html) 
  + [`r`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/r.html) 
  + [`uv_damp`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/uv_damp.html) 
  + [`wvt`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/wvt.html) 
+ Equality
  + [`isequal`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxFromFile`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Write to file
  + [`writeToFile`](/classes/nonlinear-fluxes/wvnonlinearfluxqg/writetofile.html) write information about the nonlinear flux operation to file


---