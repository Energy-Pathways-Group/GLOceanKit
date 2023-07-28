---
layout: default
title: WVTransformConstantStratification
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
---

#  WVTransformConstantStratification

Wave-vortex transformation that assumes constant stratification


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVTransformConstantStratification < <a href="/classes/wvtransform/" title="WVTransform">WVTransform</a></code></pre></div></div>

## Overview
 
  To initialization an instance of the
  WVTransformConstantStratification class you must specific the
  stratification and latitude, otherwise defaults will be assumed for
  you.
  
  ```matlab
  N0 = 3*2*pi/3600;
  wvt = WVTransformConstantStratification([100e3, 100e3, 1300],[64, 64, 65], NN0=N0,latitude=30);
  ```
 
   
  


## Topics
+ Initialization
  + [`WVTransformConstantStratification`](/classes/transforms/wvtransformconstantstratification/wvtransformconstantstratification.html) initialze a wave-vortex transform with constant stratification
  + [`waveVortexTransformFromFile`](/classes/transforms/wvtransformconstantstratification/wavevortextransformfromfile.html) Initialize a WVTransformConstantStratification instance from an existing file
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformconstantstratification/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes/transforms/wvtransformconstantstratification/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes/transforms/wvtransformconstantstratification/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes/transforms/wvtransformconstantstratification/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes/transforms/wvtransformconstantstratification/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`DCT`](/classes/transforms/wvtransformconstantstratification/dct.html) 
  + [`DFT`](/classes/transforms/wvtransformconstantstratification/dft.html) 
  + [`DST`](/classes/transforms/wvtransformconstantstratification/dst.html) 
  + [`F`](/classes/transforms/wvtransformconstantstratification/f.html) 
  + [`G`](/classes/transforms/wvtransformconstantstratification/g.html) 
  + [`N0`](/classes/transforms/wvtransformconstantstratification/n0.html) 
  + [`N2`](/classes/transforms/wvtransformconstantstratification/n2.html) 
  + [`N2AtDepth`](/classes/transforms/wvtransformconstantstratification/n2atdepth.html) 
  + [`PlaceParticlesOnIsopycnal`](/classes/transforms/wvtransformconstantstratification/placeparticlesonisopycnal.html) MAS 1/10/18 - added intext ('int' or 'both') to give option of using int vs. int+ext fields for rho_prime
  + [`ProfileTransforms`](/classes/transforms/wvtransformconstantstratification/profiletransforms.html) 
  + [`RhoBarAtDepth`](/classes/transforms/wvtransformconstantstratification/rhobaratdepth.html) 
  + [`buildTransformationMatrices`](/classes/transforms/wvtransformconstantstratification/buildtransformationmatrices.html) We renormalization the transformation matrices to directly
  + [`cg_x`](/classes/transforms/wvtransformconstantstratification/cg_x.html) 
  + [`cg_y`](/classes/transforms/wvtransformconstantstratification/cg_y.html) 
  + [`cg_z`](/classes/transforms/wvtransformconstantstratification/cg_z.html) 
  + [`diffZF`](/classes/transforms/wvtransformconstantstratification/diffzf.html) 
  + [`diffZG`](/classes/transforms/wvtransformconstantstratification/diffzg.html) 
  + [`h`](/classes/transforms/wvtransformconstantstratification/h.html) all subclasses need to have a function that returns the eigendepths
  + [`iDCT`](/classes/transforms/wvtransformconstantstratification/idct.html) 
  + [`iDFT`](/classes/transforms/wvtransformconstantstratification/idft.html) 
  + [`iDST`](/classes/transforms/wvtransformconstantstratification/idst.html) 
  + [`isHydrostatic`](/classes/transforms/wvtransformconstantstratification/ishydrostatic.html) 
  + [`rhobar`](/classes/transforms/wvtransformconstantstratification/rhobar.html) 
  + [`transformFromSpatialDomainWithF`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithf.html) 
  + [`transformFromSpatialDomainWithF_FFT`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithf_fft.html) 
  + [`transformFromSpatialDomainWithF_MM`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithf_mm.html) 
  + [`transformFromSpatialDomainWithG`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithg.html) 
  + [`transformFromSpatialDomainWithG_FFT`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithg_fft.html) df = 1/(2*(Nz-1)*dz)
  + [`transformFromSpatialDomainWithG_MM`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithg_mm.html) df = 1/(2*(Nz-1)*dz)
  + [`transformToSpatialDomainWithF`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithFAllDerivatives_FFT`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives_fft.html) 
  + [`transformToSpatialDomainWithFAllDerivatives_MM`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives_mm.html) 
  + [`transformToSpatialDomainWithF_FFT`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithf_fft.html) 
  + [`transformToSpatialDomainWithF_MM`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithf_mm.html) All coefficients are subsumbed into the transform
  + [`transformToSpatialDomainWithG`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives.html) 
  + [`transformToSpatialDomainWithGAllDerivatives_FFT`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives_fft.html) 
  + [`transformToSpatialDomainWithGAllDerivatives_MM`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives_mm.html) 
  + [`transformToSpatialDomainWithG_FFT`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithg_fft.html) 
  + [`transformToSpatialDomainWithG_MM`](/classes/transforms/wvtransformconstantstratification/transformtospatialdomainwithg_mm.html) All coefficients are subsumbed into the transform
  + [`uMaxGNormRatioForWave`](/classes/transforms/wvtransformconstantstratification/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
+ Internal
  + [`defaultPropertyAnnotations`](/classes/transforms/wvtransformconstantstratification/defaultpropertyannotations.html) return array of WVPropertyAnnotation initialized by default
+ Write to file
  + [`writeToFile`](/classes/transforms/wvtransformconstantstratification/writetofile.html) Output the `WVTransformConstantStratification` instance to file.


---