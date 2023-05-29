---
layout: default
title: WVTransformHydrostatic
parent: WV transform & model
has_children: false
has_toc: false
mathjax: true
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
  + [`WVTransformHydrostatic`](/classes-transform-and-model/wvtransformhydrostatic/wvtransformhydrostatic.html) create a wave-vortex transform for variable stratification
  + [`waveVortexTransformWithResolution`](/classes-transform-and-model/wvtransformhydrostatic/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Other
  + [`A0_HKE_factor`](/classes-transform-and-model/wvtransformhydrostatic/a0_hke_factor.html) 
  + [`A0_PE_factor`](/classes-transform-and-model/wvtransformhydrostatic/a0_pe_factor.html) 
  + [`A0_TE_factor`](/classes-transform-and-model/wvtransformhydrostatic/a0_te_factor.html) 
  + [`Apm_TE_factor`](/classes-transform-and-model/wvtransformhydrostatic/apm_te_factor.html) These convert the coefficients to their depth integrated energies
  + [`BuildProjectionOperators`](/classes-transform-and-model/wvtransformhydrostatic/buildprojectionoperators.html) Now go compute the appropriate number of modes at the
  + [`N2`](/classes-transform-and-model/wvtransformhydrostatic/n2.html) 
  + [`N2Function`](/classes-transform-and-model/wvtransformhydrostatic/n2function.html) 
  + [`P`](/classes-transform-and-model/wvtransformhydrostatic/p.html) Preconditioner for F, size(P)=[1 1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
  + [`PF`](/classes-transform-and-model/wvtransformhydrostatic/pf.html) size(PF,PG)=[Nj x Nz]
  + [`PFinv`](/classes-transform-and-model/wvtransformhydrostatic/pfinv.html) Transformation matrices
  + [`PFinvInterp`](/classes-transform-and-model/wvtransformhydrostatic/pfinvinterp.html) 
  + [`Q`](/classes-transform-and-model/wvtransformhydrostatic/q.html) Preconditioner for G, size(Q)=[1 1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
  + [`QG`](/classes-transform-and-model/wvtransformhydrostatic/qg.html) 
  + [`QGinv`](/classes-transform-and-model/wvtransformhydrostatic/qginv.html) 
  + [`QGinvInterp`](/classes-transform-and-model/wvtransformhydrostatic/qginvinterp.html) 
  + [`SetProjectionOperators`](/classes-transform-and-model/wvtransformhydrostatic/setprojectionoperators.html) 
  + [`buildInterpolationProjectionOperators`](/classes-transform-and-model/wvtransformhydrostatic/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes-transform-and-model/wvtransformhydrostatic/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`dLnN2`](/classes-transform-and-model/wvtransformhydrostatic/dlnn2.html) 
  + [`dLnN2Function`](/classes-transform-and-model/wvtransformhydrostatic/dlnn2function.html) 
  + [`diffZF`](/classes-transform-and-model/wvtransformhydrostatic/diffzf.html) 
  + [`h`](/classes-transform-and-model/wvtransformhydrostatic/h.html) [1 1 Nj]
  + [`internalModes`](/classes-transform-and-model/wvtransformhydrostatic/internalmodes.html) 
  + [`rhoFunction`](/classes-transform-and-model/wvtransformhydrostatic/rhofunction.html) function handles
  + [`rhobar`](/classes-transform-and-model/wvtransformhydrostatic/rhobar.html) on the z-grid, size(N2) = [length(z) 1];
  + [`transformFromSpatialDomainWithF`](/classes-transform-and-model/wvtransformhydrostatic/transformfromspatialdomainwithf.html) hydrostatic modes commute with the DFT
  + [`transformFromSpatialDomainWithG`](/classes-transform-and-model/wvtransformhydrostatic/transformfromspatialdomainwithg.html) hydrostatic modes commute with the DFT
  + [`transformToSpatialDomainWithF`](/classes-transform-and-model/wvtransformhydrostatic/transformtospatialdomainwithf.html) 
  + [`transformToSpatialDomainWithFAllDerivatives`](/classes-transform-and-model/wvtransformhydrostatic/transformtospatialdomainwithfallderivatives.html) 
  + [`transformToSpatialDomainWithFInterp`](/classes-transform-and-model/wvtransformhydrostatic/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithG`](/classes-transform-and-model/wvtransformhydrostatic/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGAllDerivatives`](/classes-transform-and-model/wvtransformhydrostatic/transformtospatialdomainwithgallderivatives.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes-transform-and-model/wvtransformhydrostatic/transformtospatialdomainwithginterp.html) 
  + [`uMaxGNormRatioForWave`](/classes-transform-and-model/wvtransformhydrostatic/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`zInterp`](/classes-transform-and-model/wvtransformhydrostatic/zinterp.html) 
+ Operations
  + Transformations
    + [`FMatrix`](/classes-transform-and-model/wvtransformhydrostatic/fmatrix.html) transformation matrix $$F$$
    + [`FinvMatrix`](/classes-transform-and-model/wvtransformhydrostatic/finvmatrix.html) transformation matrix $$F^{-1}$$
+ Nonlinear flux and energy transfers
  + [`nonlinearFlux`](/classes-transform-and-model/wvtransformhydrostatic/nonlinearflux.html) returns the flux of each coefficient as determined by the nonlinear flux operation
+ Initialization (Static)
  + [`waveVortexTransformFromFile`](/classes-transform-and-model/wvtransformhydrostatic/wavevortextransformfromfile.html) Initialize a WVTransformHydrostatic instance from an existing file
+ Write to file
  + [`writeToFile`](/classes-transform-and-model/wvtransformhydrostatic/writetofile.html) Output the `WVTransform` to file.


---