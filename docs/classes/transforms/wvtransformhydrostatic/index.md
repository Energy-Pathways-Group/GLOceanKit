---
layout: default
title: WVTransformHydrostatic
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
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
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformhydrostatic/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`BuildProjectionOperators`](/classes/transforms/wvtransformhydrostatic/buildprojectionoperators.html) Now go compute the appropriate number of modes at the
  + [`N2`](/classes/transforms/wvtransformhydrostatic/n2.html) 
  + [`N2Function`](/classes/transforms/wvtransformhydrostatic/n2function.html) 
  + [`P`](/classes/transforms/wvtransformhydrostatic/p.html) Preconditioner for F, size(P)=[Nj 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
  + [`PF`](/classes/transforms/wvtransformhydrostatic/pf.html) size(PF,PG)=[Nj x Nz]
  + [`PFinv`](/classes/transforms/wvtransformhydrostatic/pfinv.html) Transformation matrices
  + [`PFinvInterp`](/classes/transforms/wvtransformhydrostatic/pfinvinterp.html) 
  + [`Q`](/classes/transforms/wvtransformhydrostatic/q.html) Preconditioner for G, size(Q)=[Nj 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
  + [`QG`](/classes/transforms/wvtransformhydrostatic/qg.html) 
  + [`QGinv`](/classes/transforms/wvtransformhydrostatic/qginv.html) 
  + [`QGinvInterp`](/classes/transforms/wvtransformhydrostatic/qginvinterp.html) 
  + [`buildInterpolationProjectionOperators`](/classes/transforms/wvtransformhydrostatic/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes/transforms/wvtransformhydrostatic/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`dLnN2`](/classes/transforms/wvtransformhydrostatic/dlnn2.html) 
  + [`dLnN2Function`](/classes/transforms/wvtransformhydrostatic/dlnn2function.html) 
  + [`dftBuffer`](/classes/transforms/wvtransformhydrostatic/dftbuffer.html) 
  + [`dftConjugateIndex`](/classes/transforms/wvtransformhydrostatic/dftconjugateindex.html) 
  + [`dftPrimaryIndex`](/classes/transforms/wvtransformhydrostatic/dftprimaryindex.html) 
  + [`diffZF`](/classes/transforms/wvtransformhydrostatic/diffzf.html) 
  + [`diffZG`](/classes/transforms/wvtransformhydrostatic/diffzg.html) 
  + [`h`](/classes/transforms/wvtransformhydrostatic/h.html) [Nj 1]
  + [`h_0`](/classes/transforms/wvtransformhydrostatic/h_0.html) [Nj 1]
  + [`h_pm`](/classes/transforms/wvtransformhydrostatic/h_pm.html) [Nj 1]
  + [`iOmega`](/classes/transforms/wvtransformhydrostatic/iomega.html) 
  + [`internalModes`](/classes/transforms/wvtransformhydrostatic/internalmodes.html) 
  + [`isHydrostatic`](/classes/transforms/wvtransformhydrostatic/ishydrostatic.html) 
  + [`isequal`](/classes/transforms/wvtransformhydrostatic/isequal.html) 
  + [`rhoFunction`](/classes/transforms/wvtransformhydrostatic/rhofunction.html) function handles
  + [`rhobar`](/classes/transforms/wvtransformhydrostatic/rhobar.html) on the z-grid, size(N2) = [length(z) 1];
  + [`transformFromSpatialDomainWithFg`](/classes/transforms/wvtransformhydrostatic/transformfromspatialdomainwithfg.html) 
  + [`transformFromSpatialDomainWithFio`](/classes/transforms/wvtransformhydrostatic/transformfromspatialdomainwithfio.html) Required for transformUVEtaToWaveVortex
  + [`transformFromSpatialDomainWithGg`](/classes/transforms/wvtransformhydrostatic/transformfromspatialdomainwithgg.html) 
  + [`transformToSpatialDomainWithF`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithf.html) Perform the vertical mode matrix multiplication
  + [`transformToSpatialDomainWithFInterp`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithG`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithg.html) Perform the vertical mode matrix multiplication
  + [`transformToSpatialDomainWithGInterp`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithginterp.html) 
  + [`transformWithG_wg`](/classes/transforms/wvtransformhydrostatic/transformwithg_wg.html) 
  + [`uMaxA0`](/classes/transforms/wvtransformhydrostatic/umaxa0.html) uMax for a geostrophic mode is uMax =(g/f)*Kh*max(F_j)*abs(A0)
  + [`uMaxGNormRatioForWave`](/classes/transforms/wvtransformhydrostatic/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`wvBuffer`](/classes/transforms/wvtransformhydrostatic/wvbuffer.html) 
  + [`wvConjugateIndex`](/classes/transforms/wvtransformhydrostatic/wvconjugateindex.html) 
  + [`zInterp`](/classes/transforms/wvtransformhydrostatic/zinterp.html) 
+ Operations
  + Transformations
    + [`FMatrix`](/classes/transforms/wvtransformhydrostatic/fmatrix.html) transformation matrix $$F$$
    + [`FinvMatrix`](/classes/transforms/wvtransformhydrostatic/finvmatrix.html) transformation matrix $$F^{-1}$$
    + [`GMatrix`](/classes/transforms/wvtransformhydrostatic/gmatrix.html) transformation matrix $$G$$
    + [`GinvMatrix`](/classes/transforms/wvtransformhydrostatic/ginvmatrix.html) transformation matrix $$G^{-1}$$
+ Initialization (Static)
  + [`waveVortexTransformFromFile`](/classes/transforms/wvtransformhydrostatic/wavevortextransformfromfile.html) Initialize a WVTransformHydrostatic instance from an existing file
+ Write to file
  + [`writeToFile`](/classes/transforms/wvtransformhydrostatic/writetofile.html) Output the `WVTransform` to file.


---