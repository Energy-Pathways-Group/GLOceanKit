---
layout: default
title: WVTransformSingleMode
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
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
  + [`WVTransformSingleMode`](/classes/transforms/wvtransformsinglemode/wvtransformsinglemode.html) create a single mode wave-vortex transform
  + [`waveVortexTransformFromFile`](/classes/transforms/wvtransformsinglemode/wavevortextransformfromfile.html) Initialize a WVTransformSingleMode instance from an existing file
  + [`waveVortexTransformWithDoubleResolution`](/classes/transforms/wvtransformsinglemode/wavevortextransformwithdoubleresolution.html) create a new WVTransform with double resolution
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformsinglemode/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes/transforms/wvtransformsinglemode/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes/transforms/wvtransformsinglemode/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes/transforms/wvtransformsinglemode/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes/transforms/wvtransformsinglemode/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`buildTransformationMatrices`](/classes/transforms/wvtransformsinglemode/buildtransformationmatrices.html) Build wavenumbers
  + [`energyFluxWithMasks`](/classes/transforms/wvtransformsinglemode/energyfluxwithmasks.html) 
  + [`enstrophyFlux`](/classes/transforms/wvtransformsinglemode/enstrophyflux.html) 
  + [`h`](/classes/transforms/wvtransformsinglemode/h.html) [1 x 1]
  + [`nonlinearFluxWithMasks`](/classes/transforms/wvtransformsinglemode/nonlinearfluxwithmasks.html) 
  + [`qgpvFlux`](/classes/transforms/wvtransformsinglemode/qgpvflux.html) 
  + [`setSSH`](/classes/transforms/wvtransformsinglemode/setssh.html) 
  + [`transformFromSpatialDomainWithF`](/classes/transforms/wvtransformsinglemode/transformfromspatialdomainwithf.html) 
  + [`transformFromSpatialDomainWithG`](/classes/transforms/wvtransformsinglemode/transformfromspatialdomainwithg.html) 
  + [`transformToSpatialDomainWithF`](/classes/transforms/wvtransformsinglemode/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes/transforms/wvtransformsinglemode/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithG`](/classes/transforms/wvtransformsinglemode/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes/transforms/wvtransformsinglemode/transformtospatialdomainwithgallderivatives.html) 
  + [`uMaxGNormRatioForWave`](/classes/transforms/wvtransformsinglemode/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
+ Write to file
  + [`writeToFile`](/classes/transforms/wvtransformsinglemode/writetofile.html) Output the `WVTransformSingleMode` instance to file.


---