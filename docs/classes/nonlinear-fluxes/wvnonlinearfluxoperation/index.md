---
layout: default
title: WVNonlinearFluxOperation
has_children: false
has_toc: false
mathjax: true
parent: Nonlinear fluxes
grand_parent: Class documentation
nav_order: 1
---

#  WVNonlinearFluxOperation

Computes the nonlinear flux for a WVTransform


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVNonlinearFluxOperation < <a href="/classes/wvoperation/" title="WVOperation">WVOperation</a></code></pre></div></div>

## Overview
  
  A WVNonlinearFluxOperation defines how energy moves between the
  wave-vortex coefficients---these are the nonlinear terms in the equations
  of motion, transformed into wave-vortex space. The most basic
  implementation is a freely evolving, unforced (and undamped) flux.
  Subclasses of WVNonlinearFluxOperation can implement custom forcing and
  damping.
  
  A WVTransform is always initialized with a default nonlinear flux
  operation, but it can be overridden with,
  
  ```matlab
  wvt.nonlinearFluxOperation = SomeNewNonlinearFluxOperation();
  ```
  
  It is very likely you will want to use a custom nonlinear flux operation
  when integrating a model. In that case you would call,
  
  ```matlab
  model = WVModel(wvt,nonlinearFlux=SomeNewNonlinearFluxOperation())
  ```
  
  When creating a subclass of WVNonlinearFluxOperation, there are several
  important notes:
  
  + The output variables *must* be at least one of {Fp,Fm,F0}, in that
  order. The properties `doesFluxAp` etc. should be appropriately set to
  match the output.
  + You may also optionally output additional variables
  that are computed as a by product of your flux calculation. Those
  variables will then be cached, and will not have to be recomputed when
  needed.
  
  


## Topics
+ Initialization
  + [`WVNonlinearFluxOperation`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/wvnonlinearfluxoperation.html) create a new nonlinear flux operation
  + [`nonlinearFluxFromFile`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
  + [`nonlinearFluxWithDoubleResolution`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
+ Properties
  + [`doesFluxA0`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/doesfluxa0.html) boolean indicating whether or not this operation returns F0
  + [`doesFluxAm`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/doesfluxam.html) boolean indicating whether or not this operation returns Fm
  + [`doesFluxAp`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/doesfluxap.html) boolean indicating whether or not this operation returns Fp
+ Equality
  + [`isequal`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/isequal.html) check for equality with another nonlinear flux operation
+ Write to file
  + [`writeToFile`](/classes/nonlinear-fluxes/wvnonlinearfluxoperation/writetofile.html) write information about the nonlinear flux operation to file


---