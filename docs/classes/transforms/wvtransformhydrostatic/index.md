---
layout: default
title: WVTransformHydrostatic
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
nav_order: 2
---

#  WVTransformHydrostatic

A class for disentangling hydrostatic waves and vortices in variable stratification


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVTransformHydrostatic < <a href="/classes/wvtransform/" title="WVTransform">WVTransform</a></code></pre></div></div>

## Overview
 
  To initialization an instance of the WVTransformHydrostatic class you
  must specific the domain size, the number of grid points and *either*
  the density profile or the stratification profile.
  
  ```matlab
  N0 = 3*2*pi/3600;
  L_gm = 1300;
  N2 = @(z) N0*N0*exp(2*z/L_gm);
  wvt = WVTransformHydrostatic([100e3, 100e3, 4000],[64, 64, 65], N2=N2,latitude=30);
  ```
 
   
  


## Topics
+ Initialization
  + [`WVTransformHydrostatic`](/classes/transforms/wvtransformhydrostatic/wvtransformhydrostatic.html) create a wave-vortex transform for variable stratification
+ Operations
  + Transformations
    + [`FMatrix`](/classes/transforms/wvtransformhydrostatic/fmatrix.html) transformation matrix $$F$$
    + [`FinvMatrix`](/classes/transforms/wvtransformhydrostatic/finvmatrix.html) transformation matrix $$F^{-1}$$
    + [`GMatrix`](/classes/transforms/wvtransformhydrostatic/gmatrix.html) transformation matrix $$G$$
    + [`GinvMatrix`](/classes/transforms/wvtransformhydrostatic/ginvmatrix.html) transformation matrix $$G^{-1}$$
+ Other
  + [`BuildProjectionOperators`](/classes/transforms/wvtransformhydrostatic/buildprojectionoperators.html) Now go compute the appropriate number of modes at the
  + [`N2`](/classes/transforms/wvtransformhydrostatic/n2.html) 
  + [`N2Function`](/classes/transforms/wvtransformhydrostatic/n2function.html) 
  + [`buildInterpolationProjectionOperators`](/classes/transforms/wvtransformhydrostatic/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes/transforms/wvtransformhydrostatic/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`dLnN2`](/classes/transforms/wvtransformhydrostatic/dlnn2.html) 
  + [`dLnN2Function`](/classes/transforms/wvtransformhydrostatic/dlnn2function.html) 
  + [`rhoFunction`](/classes/transforms/wvtransformhydrostatic/rhofunction.html) function handles
  + [`rhobar`](/classes/transforms/wvtransformhydrostatic/rhobar.html) on the z-grid, size(N2) = [length(z) 1];
  + [`transformToSpatialDomainWithFInterp`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithginterp.html) 
  + [`verticalModes`](/classes/transforms/wvtransformhydrostatic/verticalmodes.html) 
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformhydrostatic/wavevortextransformwithresolution.html) 


---