---
layout: default
title: WVTransformConstantStratification
parent: WV transform & model
has_children: false
has_toc: false
mathjax: true
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
  + [`WVTransformConstantStratification`](/classes-transform-and-model/wvtransformconstantstratification/wvtransformconstantstratification.html) initialze a wave-vortex transform with constant stratification
  + [`waveVortexTransformFromFile`](/classes-transform-and-model/wvtransformconstantstratification/wavevortextransformfromfile.html) Initialize a WVTransformConstantStratification instance from an existing file
  + [`waveVortexTransformWithResolution`](/classes-transform-and-model/wvtransformconstantstratification/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes-transform-and-model/wvtransformconstantstratification/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes-transform-and-model/wvtransformconstantstratification/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes-transform-and-model/wvtransformconstantstratification/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes-transform-and-model/wvtransformconstantstratification/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`DCT`](/classes-transform-and-model/wvtransformconstantstratification/dct.html) 
  + [`DFT`](/classes-transform-and-model/wvtransformconstantstratification/dft.html) 
  + [`DST`](/classes-transform-and-model/wvtransformconstantstratification/dst.html) 
  + [`F`](/classes-transform-and-model/wvtransformconstantstratification/f.html) 
  + [`G`](/classes-transform-and-model/wvtransformconstantstratification/g.html) 
  + [`N0`](/classes-transform-and-model/wvtransformconstantstratification/n0.html) 
  + [`N2`](/classes-transform-and-model/wvtransformconstantstratification/n2.html) 
  + [`N2AtDepth`](/classes-transform-and-model/wvtransformconstantstratification/n2atdepth.html) 
  + [`PlaceParticlesOnIsopycnal`](/classes-transform-and-model/wvtransformconstantstratification/placeparticlesonisopycnal.html) MAS 1/10/18 - added intext ('int' or 'both') to give option of using int vs. int+ext fields for rho_prime
  + [`ProfileTransforms`](/classes-transform-and-model/wvtransformconstantstratification/profiletransforms.html) 
  + [`RhoBarAtDepth`](/classes-transform-and-model/wvtransformconstantstratification/rhobaratdepth.html) 
  + [`buildTransformationMatrices`](/classes-transform-and-model/wvtransformconstantstratification/buildtransformationmatrices.html) We renormalization the transformation matrices to directly
  + [`cg_x`](/classes-transform-and-model/wvtransformconstantstratification/cg_x.html) 
  + [`cg_y`](/classes-transform-and-model/wvtransformconstantstratification/cg_y.html) 
  + [`cg_z`](/classes-transform-and-model/wvtransformconstantstratification/cg_z.html) 
  + [`diffZF`](/classes-transform-and-model/wvtransformconstantstratification/diffzf.html) 
  + [`diffZG`](/classes-transform-and-model/wvtransformconstantstratification/diffzg.html) 
  + [`h`](/classes-transform-and-model/wvtransformconstantstratification/h.html) all subclasses need to have a function that returns the eigendepths
  + [`iDCT`](/classes-transform-and-model/wvtransformconstantstratification/idct.html) 
  + [`iDFT`](/classes-transform-and-model/wvtransformconstantstratification/idft.html) 
  + [`iDST`](/classes-transform-and-model/wvtransformconstantstratification/idst.html) 
  + [`isHydrostatic`](/classes-transform-and-model/wvtransformconstantstratification/ishydrostatic.html) 
  + [`rhobar`](/classes-transform-and-model/wvtransformconstantstratification/rhobar.html) 
  + [`transformFromSpatialDomainWithF`](/classes-transform-and-model/wvtransformconstantstratification/transformfromspatialdomainwithf.html) 
  + [`transformFromSpatialDomainWithF_FFT`](/classes-transform-and-model/wvtransformconstantstratification/transformfromspatialdomainwithf_fft.html) 
  + [`transformFromSpatialDomainWithF_MM`](/classes-transform-and-model/wvtransformconstantstratification/transformfromspatialdomainwithf_mm.html) 
  + [`transformFromSpatialDomainWithG`](/classes-transform-and-model/wvtransformconstantstratification/transformfromspatialdomainwithg.html) 
  + [`transformFromSpatialDomainWithG_FFT`](/classes-transform-and-model/wvtransformconstantstratification/transformfromspatialdomainwithg_fft.html) df = 1/(2*(Nz-1)*dz)
  + [`transformFromSpatialDomainWithG_MM`](/classes-transform-and-model/wvtransformconstantstratification/transformfromspatialdomainwithg_mm.html) df = 1/(2*(Nz-1)*dz)
  + [`transformToSpatialDomainWithF`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithFAllDerivatives_FFT`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives_fft.html) 
  + [`transformToSpatialDomainWithFAllDerivatives_MM`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithfallderivatives_mm.html) 
  + [`transformToSpatialDomainWithF_FFT`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithf_fft.html) 
  + [`transformToSpatialDomainWithF_MM`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithf_mm.html) All coefficients are subsumbed into the transform
  + [`transformToSpatialDomainWithG`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives.html) 
  + [`transformToSpatialDomainWithGAllDerivatives_FFT`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives_fft.html) 
  + [`transformToSpatialDomainWithGAllDerivatives_MM`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithgallderivatives_mm.html) 
  + [`transformToSpatialDomainWithG_FFT`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithg_fft.html) 
  + [`transformToSpatialDomainWithG_MM`](/classes-transform-and-model/wvtransformconstantstratification/transformtospatialdomainwithg_mm.html) All coefficients are subsumbed into the transform
  + [`uMaxGNormRatioForWave`](/classes-transform-and-model/wvtransformconstantstratification/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
+ Internal
  + [`defaultPropertyAnnotations`](/classes-transform-and-model/wvtransformconstantstratification/defaultpropertyannotations.html) return array of WVPropertyAnnotation initialized by default
+ Write to file
  + [`writeToFile`](/classes-transform-and-model/wvtransformconstantstratification/writetofile.html) Output the `WVTransformConstantStratification` instance to file.


---