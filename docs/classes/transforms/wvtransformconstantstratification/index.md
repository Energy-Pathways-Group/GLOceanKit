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
  + [`WVTransformConstantStratification`](/classes/wvtransformconstantstratification/wvtransformconstantstratification.html) initialze a wave-vortex transform with constant stratification
  + [`waveVortexTransformFromFile`](/classes/wvtransformconstantstratification/wavevortextransformfromfile.html) Initialize a WVTransformConstantStratification instance from an existing file
  + [`waveVortexTransformWithResolution`](/classes/wvtransformconstantstratification/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes/wvtransformconstantstratification/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes/wvtransformconstantstratification/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes/wvtransformconstantstratification/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes/wvtransformconstantstratification/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`DCT`](/classes/wvtransformconstantstratification/dct.html) 
  + [`DFT`](/classes/wvtransformconstantstratification/dft.html) 
  + [`DST`](/classes/wvtransformconstantstratification/dst.html) 
  + [`F`](/classes/wvtransformconstantstratification/f.html) 
  + [`G`](/classes/wvtransformconstantstratification/g.html) 
  + [`N0`](/classes/wvtransformconstantstratification/n0.html) 
  + [`N2`](/classes/wvtransformconstantstratification/n2.html) 
  + [`N2AtDepth`](/classes/wvtransformconstantstratification/n2atdepth.html) 
  + [`PlaceParticlesOnIsopycnal`](/classes/wvtransformconstantstratification/placeparticlesonisopycnal.html) MAS 1/10/18 - added intext ('int' or 'both') to give option of using int vs. int+ext fields for rho_prime
  + [`ProfileTransforms`](/classes/wvtransformconstantstratification/profiletransforms.html) 
  + [`RhoBarAtDepth`](/classes/wvtransformconstantstratification/rhobaratdepth.html) 
  + [`buildTransformationMatrices`](/classes/wvtransformconstantstratification/buildtransformationmatrices.html) We renormalization the transformation matrices to directly
  + [`cg_x`](/classes/wvtransformconstantstratification/cg_x.html) 
  + [`cg_y`](/classes/wvtransformconstantstratification/cg_y.html) 
  + [`cg_z`](/classes/wvtransformconstantstratification/cg_z.html) 
  + [`diffZF`](/classes/wvtransformconstantstratification/diffzf.html) 
  + [`diffZG`](/classes/wvtransformconstantstratification/diffzg.html) 
  + [`h`](/classes/wvtransformconstantstratification/h.html) all subclasses need to have a function that returns the eigendepths
  + [`iDCT`](/classes/wvtransformconstantstratification/idct.html) 
  + [`iDFT`](/classes/wvtransformconstantstratification/idft.html) 
  + [`iDST`](/classes/wvtransformconstantstratification/idst.html) 
  + [`isHydrostatic`](/classes/wvtransformconstantstratification/ishydrostatic.html) 
  + [`rhobar`](/classes/wvtransformconstantstratification/rhobar.html) 
  + [`transformFromSpatialDomainWithF`](/classes/wvtransformconstantstratification/transformfromspatialdomainwithf.html) 
  + [`transformFromSpatialDomainWithF_FFT`](/classes/wvtransformconstantstratification/transformfromspatialdomainwithf_fft.html) 
  + [`transformFromSpatialDomainWithF_MM`](/classes/wvtransformconstantstratification/transformfromspatialdomainwithf_mm.html) 
  + [`transformFromSpatialDomainWithG`](/classes/wvtransformconstantstratification/transformfromspatialdomainwithg.html) 
  + [`transformFromSpatialDomainWithG_FFT`](/classes/wvtransformconstantstratification/transformfromspatialdomainwithg_fft.html) df = 1/(2*(Nz-1)*dz)
  + [`transformFromSpatialDomainWithG_MM`](/classes/wvtransformconstantstratification/transformfromspatialdomainwithg_mm.html) df = 1/(2*(Nz-1)*dz)
  + [`transformToSpatialDomainWithF`](/classes/wvtransformconstantstratification/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithFAllDerivatives_FFT`](/classes/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives_fft.html) 
  + [`transformToSpatialDomainWithFAllDerivatives_MM`](/classes/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives_mm.html) 
  + [`transformToSpatialDomainWithF_FFT`](/classes/wvtransformconstantstratification/transformtospatialdomainwithf_fft.html) 
  + [`transformToSpatialDomainWithF_MM`](/classes/wvtransformconstantstratification/transformtospatialdomainwithf_mm.html) All coefficients are subsumbed into the transform
  + [`transformToSpatialDomainWithG`](/classes/wvtransformconstantstratification/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives.html) 
  + [`transformToSpatialDomainWithGAllDerivatives_FFT`](/classes/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives_fft.html) 
  + [`transformToSpatialDomainWithGAllDerivatives_MM`](/classes/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives_mm.html) 
  + [`transformToSpatialDomainWithG_FFT`](/classes/wvtransformconstantstratification/transformtospatialdomainwithg_fft.html) 
  + [`transformToSpatialDomainWithG_MM`](/classes/wvtransformconstantstratification/transformtospatialdomainwithg_mm.html) All coefficients are subsumbed into the transform
  + [`uMaxGNormRatioForWave`](/classes/wvtransformconstantstratification/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
+ Internal
  + [`defaultPropertyAnnotations`](/classes/wvtransformconstantstratification/defaultpropertyannotations.html) return array of WVPropertyAnnotation initialized by default
+ Write to file
  + [`writeToFile`](/classes/wvtransformconstantstratification/writetofile.html) Output the `WVTransformConstantStratification` instance to file.


---