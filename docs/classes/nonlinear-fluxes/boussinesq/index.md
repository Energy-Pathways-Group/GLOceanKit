---
layout: default
title: Boussinesq
has_children: false
has_toc: false
mathjax: true
parent: Nonlinear fluxes
grand_parent: Class documentation
nav_order: 2
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
  + [`AA`](/classes/boussinesq/aa.html) 
  + [`Boussinesq`](/classes/boussinesq/boussinesq.html) initialize the Boussinesq nonlinear flux
  + [`dLnN2`](/classes/boussinesq/dlnn2.html) 
  + [`damp`](/classes/boussinesq/damp.html) 
  + [`nu_xy`](/classes/boussinesq/nu_xy.html) 
  + [`nu_z`](/classes/boussinesq/nu_z.html) 
  + [`shouldAntialias`](/classes/boussinesq/shouldantialias.html) 
+ Computation
  + [`compute`](/classes/boussinesq/compute.html) the promised variable
+ Equality
  + [`isequal`](/classes/boussinesq/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxFromFile`](/classes/boussinesq/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes/boussinesq/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Write to file
  + [`writeToFile`](/classes/boussinesq/writetofile.html) write information about the nonlinear flux operation to file


---