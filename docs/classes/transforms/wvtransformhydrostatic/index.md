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
  + [`WVTransformHydrostatic`](/classes/wvtransformhydrostatic/wvtransformhydrostatic.html) create a wave-vortex transform for variable stratification
  + [`waveVortexTransformWithResolution`](/classes/wvtransformhydrostatic/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes/wvtransformhydrostatic/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes/wvtransformhydrostatic/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes/wvtransformhydrostatic/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes/wvtransformhydrostatic/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`BuildProjectionOperators`](/classes/wvtransformhydrostatic/buildprojectionoperators.html) Now go compute the appropriate number of modes at the
  + [`N2`](/classes/wvtransformhydrostatic/n2.html) 
  + [`N2Function`](/classes/wvtransformhydrostatic/n2function.html) 
  + [`P`](/classes/wvtransformhydrostatic/p.html) Preconditioner for F, size(P)=[1 1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
  + [`PF`](/classes/wvtransformhydrostatic/pf.html) size(PF,PG)=[Nj x Nz]
  + [`PFinv`](/classes/wvtransformhydrostatic/pfinv.html) Transformation matrices
  + [`PFinvInterp`](/classes/wvtransformhydrostatic/pfinvinterp.html) 
  + [`Q`](/classes/wvtransformhydrostatic/q.html) Preconditioner for G, size(Q)=[1 1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
  + [`QG`](/classes/wvtransformhydrostatic/qg.html) 
  + [`QGinv`](/classes/wvtransformhydrostatic/qginv.html) 
  + [`QGinvInterp`](/classes/wvtransformhydrostatic/qginvinterp.html) 
  + [`SetProjectionOperators`](/classes/wvtransformhydrostatic/setprojectionoperators.html) 
  + [`buildInterpolationProjectionOperators`](/classes/wvtransformhydrostatic/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes/wvtransformhydrostatic/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`dLnN2`](/classes/wvtransformhydrostatic/dlnn2.html) 
  + [`dLnN2Function`](/classes/wvtransformhydrostatic/dlnn2function.html) 
  + [`diffZF`](/classes/wvtransformhydrostatic/diffzf.html) 
  + [`diffZG`](/classes/wvtransformhydrostatic/diffzg.html) 
  + [`h`](/classes/wvtransformhydrostatic/h.html) [1 1 Nj]
  + [`internalModes`](/classes/wvtransformhydrostatic/internalmodes.html) 
  + [`rhoFunction`](/classes/wvtransformhydrostatic/rhofunction.html) function handles
  + [`rhobar`](/classes/wvtransformhydrostatic/rhobar.html) on the z-grid, size(N2) = [length(z) 1];
  + [`transformFromSpatialDomainWithF`](/classes/wvtransformhydrostatic/transformfromspatialdomainwithf.html) hydrostatic modes commute with the DFT
  + [`transformFromSpatialDomainWithG`](/classes/wvtransformhydrostatic/transformfromspatialdomainwithg.html) hydrostatic modes commute with the DFT
  + [`transformToSpatialDomainWithF`](/classes/wvtransformhydrostatic/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wvtransformhydrostatic/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithFInterp`](/classes/wvtransformhydrostatic/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithG`](/classes/wvtransformhydrostatic/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wvtransformhydrostatic/transformtospatialdomainwithgallderivatives.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes/wvtransformhydrostatic/transformtospatialdomainwithginterp.html) 
  + [`uMaxGNormRatioForWave`](/classes/wvtransformhydrostatic/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`zInterp`](/classes/wvtransformhydrostatic/zinterp.html) 
+ Operations
  + Transformations
    + [`FMatrix`](/classes/wvtransformhydrostatic/fmatrix.html) transformation matrix $$F$$
    + [`FinvMatrix`](/classes/wvtransformhydrostatic/finvmatrix.html) transformation matrix $$F^{-1}$$
+ Nonlinear flux and energy transfers
  + [`nonlinearFlux`](/classes/wvtransformhydrostatic/nonlinearflux.html) returns the flux of each coefficient as determined by the nonlinear flux operation
+ Initialization (Static)
  + [`waveVortexTransformFromFile`](/classes/wvtransformhydrostatic/wavevortextransformfromfile.html) Initialize a WVTransformHydrostatic instance from an existing file
+ Write to file
  + [`writeToFile`](/classes/wvtransformhydrostatic/writetofile.html) Output the `WVTransform` to file.


---