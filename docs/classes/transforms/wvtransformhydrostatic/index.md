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
  + [`A0_HKE_factor`](/classes/transforms/wvtransformhydrostatic/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes/transforms/wvtransformhydrostatic/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes/transforms/wvtransformhydrostatic/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes/transforms/wvtransformhydrostatic/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`BuildProjectionOperators`](/classes/transforms/wvtransformhydrostatic/buildprojectionoperators.html) Now go compute the appropriate number of modes at the
  + [`N2`](/classes/transforms/wvtransformhydrostatic/n2.html) 
  + [`N2Function`](/classes/transforms/wvtransformhydrostatic/n2function.html) 
  + [`P`](/classes/transforms/wvtransformhydrostatic/p.html) Preconditioner for F, size(P)=[1 1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
  + [`PF`](/classes/transforms/wvtransformhydrostatic/pf.html) size(PF,PG)=[Nj x Nz]
  + [`PFinv`](/classes/transforms/wvtransformhydrostatic/pfinv.html) Transformation matrices
  + [`PFinvInterp`](/classes/transforms/wvtransformhydrostatic/pfinvinterp.html) 
  + [`Q`](/classes/transforms/wvtransformhydrostatic/q.html) Preconditioner for G, size(Q)=[1 1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
  + [`QG`](/classes/transforms/wvtransformhydrostatic/qg.html) 
  + [`QGinv`](/classes/transforms/wvtransformhydrostatic/qginv.html) 
  + [`QGinvInterp`](/classes/transforms/wvtransformhydrostatic/qginvinterp.html) 
  + [`SetProjectionOperators`](/classes/transforms/wvtransformhydrostatic/setprojectionoperators.html) 
  + [`buildInterpolationProjectionOperators`](/classes/transforms/wvtransformhydrostatic/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes/transforms/wvtransformhydrostatic/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`dLnN2`](/classes/transforms/wvtransformhydrostatic/dlnn2.html) 
  + [`dLnN2Function`](/classes/transforms/wvtransformhydrostatic/dlnn2function.html) 
  + [`diffZF`](/classes/transforms/wvtransformhydrostatic/diffzf.html) 
  + [`diffZG`](/classes/transforms/wvtransformhydrostatic/diffzg.html) 
  + [`h`](/classes/transforms/wvtransformhydrostatic/h.html) [1 1 Nj]
  + [`internalModes`](/classes/transforms/wvtransformhydrostatic/internalmodes.html) 
  + [`rhoFunction`](/classes/transforms/wvtransformhydrostatic/rhofunction.html) function handles
  + [`rhobar`](/classes/transforms/wvtransformhydrostatic/rhobar.html) on the z-grid, size(N2) = [length(z) 1];
  + [`transformFromSpatialDomainWithF`](/classes/transforms/wvtransformhydrostatic/transformfromspatialdomainwithf.html) hydrostatic modes commute with the DFT
  + [`transformFromSpatialDomainWithG`](/classes/transforms/wvtransformhydrostatic/transformfromspatialdomainwithg.html) hydrostatic modes commute with the DFT
  + [`transformToSpatialDomainWithF`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithFInterp`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithG`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithgallderivatives.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes/transforms/wvtransformhydrostatic/transformtospatialdomainwithginterp.html) 
  + [`uMaxGNormRatioForWave`](/classes/transforms/wvtransformhydrostatic/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`zInterp`](/classes/transforms/wvtransformhydrostatic/zinterp.html) 
+ Operations
  + Transformations
    + [`FMatrix`](/classes/transforms/wvtransformhydrostatic/fmatrix.html) transformation matrix $$F$$
    + [`FinvMatrix`](/classes/transforms/wvtransformhydrostatic/finvmatrix.html) transformation matrix $$F^{-1}$$
+ Nonlinear flux and energy transfers
  + [`nonlinearFlux`](/classes/transforms/wvtransformhydrostatic/nonlinearflux.html) returns the flux of each coefficient as determined by the nonlinear flux operation
+ Initialization (Static)
  + [`waveVortexTransformFromFile`](/classes/transforms/wvtransformhydrostatic/wavevortextransformfromfile.html) Initialize a WVTransformHydrostatic instance from an existing file
+ Write to file
  + [`writeToFile`](/classes/transforms/wvtransformhydrostatic/writetofile.html) Output the `WVTransform` to file.


---