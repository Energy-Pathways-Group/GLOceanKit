---
layout: default
title: QGPVE
has_children: false
has_toc: false
mathjax: true
parent: Nonlinear fluxes
grand_parent: Class documentation
nav_order: 4
---

#  QGPVE

3D quasigeostrophic potential vorticity flux


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>QGPVE < <a href="/classes/wvnonlinearfluxoperation/" title="WVNonlinearFluxOperation">WVNonlinearFluxOperation</a></code></pre></div></div>

## Overview
 
  The 3D quasigeostrophic potential vorticity flux will only use and
  modify the A0 coefficients.
 
  To initialize the QGPVE,
 
  ```matlab
  model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));
  ```
 
    


## Topics
+ Initializing
+ Other
  + [`A0PV`](/classes/nonlinear-fluxes/qgpve/a0pv.html) conversion from PV to A0
  + [`PVA0`](/classes/nonlinear-fluxes/qgpve/pva0.html) conversion from A0 to PV
  + [`QGPVE`](/classes/nonlinear-fluxes/qgpve/qgpve.html) initialize 3D quasigeostrophic potential vorticity flux
  + [`RVA0`](/classes/nonlinear-fluxes/qgpve/rva0.html) conversion from A0 to RV
  + [`beta`](/classes/nonlinear-fluxes/qgpve/beta.html) 
  + [`compute`](/classes/nonlinear-fluxes/qgpve/compute.html) Apply operator S---defined in (C4) in the manuscript
  + [`damp`](/classes/nonlinear-fluxes/qgpve/damp.html) 
  + [`dampingTimeScale`](/classes/nonlinear-fluxes/qgpve/dampingtimescale.html) 
  + [`nu`](/classes/nonlinear-fluxes/qgpve/nu.html) 
  + [`r`](/classes/nonlinear-fluxes/qgpve/r.html) 
  + [`uEady`](/classes/nonlinear-fluxes/qgpve/ueady.html) 
+ Equality
  + [`isequal`](/classes/nonlinear-fluxes/qgpve/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxFromFile`](/classes/nonlinear-fluxes/qgpve/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes/nonlinear-fluxes/qgpve/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Write to file
  + [`writeToFile`](/classes/nonlinear-fluxes/qgpve/writetofile.html) write information about the nonlinear flux operation to file


---