---
layout: default
title: WVTransformSingleMode
parent: WV transform & model
has_children: false
has_toc: false
mathjax: true
---

#  WVTransformSingleMode

A class for disentangling waves and vortices in a single layer


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVTransformSingleMode < <a href="/classes/wvtransform/" title="WVTransform">WVTransform</a></code></pre></div></div>

## Overview
 
  This is a two-dimensional, single-layer which may be interepreted as
  the sea-surface height. The 'h' parameter is the equivalent depth,
  and 0.80 m is a typical value for the first baroclinic mode.
 
  ```matlab
  Lxy = 50e3;
  Nxy = 256;
  latitude = 25;
  wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);
  ```
 
   
  


## Topics
+ Initialization
  + [`WVTransformSingleMode`](/classes-transform-and-model/wvtransformsinglemode/wvtransformsinglemode.html) create a single mode wave-vortex transform
  + [`waveVortexTransformFromFile`](/classes-transform-and-model/wvtransformsinglemode/wavevortextransformfromfile.html) Initialize a WVTransformSingleMode instance from an existing file
  + [`waveVortexTransformWithDoubleResolution`](/classes-transform-and-model/wvtransformsinglemode/wavevortextransformwithdoubleresolution.html) create a new WVTransform with double resolution
  + [`waveVortexTransformWithResolution`](/classes-transform-and-model/wvtransformsinglemode/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes-transform-and-model/wvtransformsinglemode/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes-transform-and-model/wvtransformsinglemode/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes-transform-and-model/wvtransformsinglemode/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes-transform-and-model/wvtransformsinglemode/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`buildTransformationMatrices`](/classes-transform-and-model/wvtransformsinglemode/buildtransformationmatrices.html) Build wavenumbers
  + [`energyFluxWithMasks`](/classes-transform-and-model/wvtransformsinglemode/energyfluxwithmasks.html) 
  + [`enstrophyFlux`](/classes-transform-and-model/wvtransformsinglemode/enstrophyflux.html) 
  + [`h`](/classes-transform-and-model/wvtransformsinglemode/h.html) [1 x 1]
  + [`nonlinearFluxWithMasks`](/classes-transform-and-model/wvtransformsinglemode/nonlinearfluxwithmasks.html) 
  + [`qgpvFlux`](/classes-transform-and-model/wvtransformsinglemode/qgpvflux.html) 
  + [`setSSH`](/classes-transform-and-model/wvtransformsinglemode/setssh.html) 
  + [`transformFromSpatialDomainWithF`](/classes-transform-and-model/wvtransformsinglemode/transformfromspatialdomainwithf.html) 
  + [`transformFromSpatialDomainWithG`](/classes-transform-and-model/wvtransformsinglemode/transformfromspatialdomainwithg.html) 
  + [`transformToSpatialDomainWithF`](/classes-transform-and-model/wvtransformsinglemode/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes-transform-and-model/wvtransformsinglemode/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithG`](/classes-transform-and-model/wvtransformsinglemode/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes-transform-and-model/wvtransformsinglemode/transformtospatialdomainwithgallderivatives.html) 
  + [`uMaxGNormRatioForWave`](/classes-transform-and-model/wvtransformsinglemode/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
+ Write to file
  + [`writeToFile`](/classes-transform-and-model/wvtransformsinglemode/writetofile.html) Output the `WVTransformSingleMode` instance to file.


---