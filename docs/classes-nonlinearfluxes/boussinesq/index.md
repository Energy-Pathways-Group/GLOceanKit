---
layout: default
title: Boussinesq
parent: Nonlinear flux operations
has_children: false
has_toc: false
mathjax: true
---

#  Boussinesq

3D nonlinear flux for Boussinesq flow, appropriate for numerical modeling


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>Boussinesq < <a href="/classes/wvnonlinearfluxoperation/" title="WVNonlinearFluxOperation">WVNonlinearFluxOperation</a></code></pre></div></div>

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
  model = WVModel(wvt,nonlinearFlux=Boussinesq(wvt,shouldAntialias=1,uv_damp=wvt.uMax));
  ```
 
    


## Topics
+ Initializing
+ Other
  + [`AA`](/classes-nonlinearfluxes/boussinesq/aa.html) 
  + [`Boussinesq`](/classes-nonlinearfluxes/boussinesq/boussinesq.html) initialize the Boussinesq nonlinear flux
  + [`dLnN2`](/classes-nonlinearfluxes/boussinesq/dlnn2.html) 
  + [`damp`](/classes-nonlinearfluxes/boussinesq/damp.html) 
  + [`nu_xy`](/classes-nonlinearfluxes/boussinesq/nu_xy.html) 
  + [`nu_z`](/classes-nonlinearfluxes/boussinesq/nu_z.html) 
  + [`shouldAntialias`](/classes-nonlinearfluxes/boussinesq/shouldantialias.html) 
+ Computation
  + [`compute`](/classes-nonlinearfluxes/boussinesq/compute.html) the promised variable
+ Equality
  + [`isequal`](/classes-nonlinearfluxes/boussinesq/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxFromFile`](/classes-nonlinearfluxes/boussinesq/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes-nonlinearfluxes/boussinesq/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Write to file
  + [`writeToFile`](/classes-nonlinearfluxes/boussinesq/writetofile.html) write information about the nonlinear flux operation to file


---