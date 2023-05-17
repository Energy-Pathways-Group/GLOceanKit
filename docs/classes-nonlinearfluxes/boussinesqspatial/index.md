---
layout: default
title: BoussinesqSpatial
parent: Nonlinear flux operations
has_children: false
has_toc: false
mathjax: true
---

#  BoussinesqSpatial

3D nonlinear flux for Boussinesq flow, computed in the spatial domain


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>BoussinesqSpatial < <a href="/classes/wvnonlinearfluxoperation/" title="WVNonlinearFluxOperation">WVNonlinearFluxOperation</a></code></pre></div></div>

## Overview
 
  Computes the nonlinear flux for a Boussinesq model. This class is not
  intended to be used for numerical modeling as it does not have any
  antialiasing or damping, but is indended as an example. The
  implementation is *simple* and follows directly from the equations of
  motion, but it is not the fastest implementation. To compute
  nonlinear fluxes appropriate for numerical modeling, use the
  [Boussinesq](/classes/boussinesq/) class.
 
    


## Topics
+ Initializing
  + [`BoussinesqSpatial`](/classes-nonlinearfluxes/boussinesqspatial/boussinesqspatial.html) 3D nonlinear flux for Boussinesq flow, computed in the spatial domain
+ Computation
  + [`compute`](/classes-nonlinearfluxes/boussinesqspatial/compute.html) the promised variable
+ Other
  + [`dLnN2`](/classes-nonlinearfluxes/boussinesqspatial/dlnn2.html) 


---