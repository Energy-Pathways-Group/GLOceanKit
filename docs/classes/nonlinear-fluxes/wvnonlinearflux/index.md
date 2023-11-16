---
layout: default
title: WVNonlinearFlux
has_children: false
has_toc: false
mathjax: true
parent: Nonlinear fluxes
grand_parent: Class documentation
nav_order: 2
---

#  WVNonlinearFlux

3D nonlinear flux for Boussinesq flow, appropriate for numerical modeling


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>WVNonlinearFlux < <a href="/classes/wvnonlinearfluxoperation/" title="WVNonlinearFluxOperation">WVNonlinearFluxOperation</a></code></pre></div></div>

## Overview
 
  Computes the nonlinear flux for a Boussinesq model, and has options
  for anti-aliasing and damping appropriate for running a numerical
  model. This is *not* the simplest implementation, but instead adds
  some complexity in favor of speed. The [BoussinesqSpatial](/classes/boussinesqspatial/) class
  shows a simple implementation.
 
  The damping is a simple Laplacian, but with a spectral vanishing
  viscosity (SVV) operator applied that prevents any damping below a
  cutoff wavenumber. The SVV operator adjusts the wavenumbers being
  damped depending on whether anti-aliasing is applied.
 
  This is most often used when initializing a model, e.g.,
 
  ```matlab
  model = WVModel(wvt,nonlinearFlux=WVNonlinearFlux(wvt,shouldAntialias=1,uv_damp=wvt.uMax));
  ```
 
    


## Topics
+ Initializing
+ Other
  + [`AA`](/classes/nonlinear-fluxes/wvnonlinearflux/aa.html) 
  + [`WVNonlinearFlux`](/classes/nonlinear-fluxes/wvnonlinearflux/wvnonlinearflux.html) initialize the WVNonlinearFlux nonlinear flux
  + [`beta`](/classes/nonlinear-fluxes/wvnonlinearflux/beta.html) 
  + [`betaA0`](/classes/nonlinear-fluxes/wvnonlinearflux/betaa0.html) 
  + [`buildDampingOperator`](/classes/nonlinear-fluxes/wvnonlinearflux/builddampingoperator.html) 
  + [`dLnN2`](/classes/nonlinear-fluxes/wvnonlinearflux/dlnn2.html) 
  + [`damp`](/classes/nonlinear-fluxes/wvnonlinearflux/damp.html) 
  + [`dampingTimeScale`](/classes/nonlinear-fluxes/wvnonlinearflux/dampingtimescale.html) 
  + [`k_damp`](/classes/nonlinear-fluxes/wvnonlinearflux/k_damp.html) wavenumber at which the first small scale damping starts.
  + [`nu_xy`](/classes/nonlinear-fluxes/wvnonlinearflux/nu_xy.html) 
  + [`nu_z`](/classes/nonlinear-fluxes/wvnonlinearflux/nu_z.html) 
  + [`r`](/classes/nonlinear-fluxes/wvnonlinearflux/r.html) 
  + [`shouldAntialias`](/classes/nonlinear-fluxes/wvnonlinearflux/shouldantialias.html) 
  + [`spatialFlux`](/classes/nonlinear-fluxes/wvnonlinearflux/spatialflux.html) a subclass can override this, and then modify the spatial
  + [`uv_damp`](/classes/nonlinear-fluxes/wvnonlinearflux/uv_damp.html) 
  + [`wvt`](/classes/nonlinear-fluxes/wvnonlinearflux/wvt.html) 
+ Computation
  + [`compute`](/classes/nonlinear-fluxes/wvnonlinearflux/compute.html) the promised variable
+ Equality
  + [`isequal`](/classes/nonlinear-fluxes/wvnonlinearflux/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxFromFile`](/classes/nonlinear-fluxes/wvnonlinearflux/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes/nonlinear-fluxes/wvnonlinearflux/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Write to file
  + [`writeToFile`](/classes/nonlinear-fluxes/wvnonlinearflux/writetofile.html) write information about the nonlinear flux operation to file


---