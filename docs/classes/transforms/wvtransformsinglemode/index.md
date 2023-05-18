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
  + [`WVTransformSingleMode`](/classes/wvtransformsinglemode/wvtransformsinglemode.html) create a single mode wave-vortex transform
  + [`waveVortexTransformFromFile`](/classes/wvtransformsinglemode/wavevortextransformfromfile.html) Initialize a WVTransformSingleMode instance from an existing file
  + [`waveVortexTransformWithDoubleResolution`](/classes/wvtransformsinglemode/wavevortextransformwithdoubleresolution.html) create a new WVTransform with double resolution
  + [`waveVortexTransformWithResolution`](/classes/wvtransformsinglemode/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes/wvtransformsinglemode/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes/wvtransformsinglemode/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes/wvtransformsinglemode/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes/wvtransformsinglemode/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`buildTransformationMatrices`](/classes/wvtransformsinglemode/buildtransformationmatrices.html) Build wavenumbers
  + [`energyFluxWithMasks`](/classes/wvtransformsinglemode/energyfluxwithmasks.html) 
  + [`enstrophyFlux`](/classes/wvtransformsinglemode/enstrophyflux.html) 
  + [`h`](/classes/wvtransformsinglemode/h.html) [1 x 1]
  + [`nonlinearFluxWithMasks`](/classes/wvtransformsinglemode/nonlinearfluxwithmasks.html) 
  + [`qgpvFlux`](/classes/wvtransformsinglemode/qgpvflux.html) 
  + [`setSSH`](/classes/wvtransformsinglemode/setssh.html) 
  + [`transformFromSpatialDomainWithF`](/classes/wvtransformsinglemode/transformfromspatialdomainwithf.html) 
  + [`transformFromSpatialDomainWithG`](/classes/wvtransformsinglemode/transformfromspatialdomainwithg.html) 
  + [`transformToSpatialDomainWithF`](/classes/wvtransformsinglemode/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wvtransformsinglemode/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithG`](/classes/wvtransformsinglemode/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wvtransformsinglemode/transformtospatialdomainwithgallderivatives.html) 
  + [`uMaxGNormRatioForWave`](/classes/wvtransformsinglemode/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
+ Write to file
  + [`writeToFile`](/classes/wvtransformsinglemode/writetofile.html) Output the `WVTransformSingleMode` instance to file.


---