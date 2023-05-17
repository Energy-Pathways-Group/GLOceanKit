---
layout: default
title: QGPVE
parent: Nonlinear flux operations
has_children: false
has_toc: false
mathjax: true
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
  + [`A0PV`](/classes-nonlinearfluxes/qgpve/a0pv.html) conversion from PV to A0
  + [`PVA0`](/classes-nonlinearfluxes/qgpve/pva0.html) conversion from A0 to PV
  + [`QGPVE`](/classes-nonlinearfluxes/qgpve/qgpve.html) initialize 3D quasigeostrophic potential vorticity flux
  + [`RVA0`](/classes-nonlinearfluxes/qgpve/rva0.html) conversion from A0 to RV
  + [`beta`](/classes-nonlinearfluxes/qgpve/beta.html) 
  + [`compute`](/classes-nonlinearfluxes/qgpve/compute.html) Apply operator S---defined in (C4) in the manuscript
  + [`damp`](/classes-nonlinearfluxes/qgpve/damp.html) 
  + [`dampingTimeScale`](/classes-nonlinearfluxes/qgpve/dampingtimescale.html) 
  + [`nu`](/classes-nonlinearfluxes/qgpve/nu.html) 
  + [`r`](/classes-nonlinearfluxes/qgpve/r.html) 
  + [`uEady`](/classes-nonlinearfluxes/qgpve/ueady.html) 
+ Equality
  + [`isequal`](/classes-nonlinearfluxes/qgpve/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxFromFile`](/classes-nonlinearfluxes/qgpve/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes-nonlinearfluxes/qgpve/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Write to file
  + [`writeToFile`](/classes-nonlinearfluxes/qgpve/writetofile.html) write information about the nonlinear flux operation to file


---