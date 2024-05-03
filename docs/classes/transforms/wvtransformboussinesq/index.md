---
layout: default
title: WVTransformBoussinesq
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
---

#  WVTransformBoussinesq

3D hydrostatic Boussinesq model with arbitrary stratification solved


---

## Overview
  in wave-vortex space
 
  Couple of different initialization paths:
  1) You want to run this as a prognostic model and therefore want
     the chebyshev points automatically found for you
        Init([Lx Ly Lz], [Nx Ny Nz], latitude, rho)
 
  2) You want to run this as a diagnostic model and therefore want
     to specify the depths and modes yourself
        Init([Lx Ly Lz], [Nx Ny Nz], latitude, rho, 'zgrid', z)


## Topics
+ Other
  + [`Aklz`](/classes/transforms/wvtransformboussinesq/aklz.html) 
  + [`BuildProjectionOperators`](/classes/transforms/wvtransformboussinesq/buildprojectionoperators.html) 
  + [`BuildProjectionOperatorsForGeostrophicModes`](/classes/transforms/wvtransformboussinesq/buildprojectionoperatorsforgeostrophicmodes.html) Now go compute the appropriate number of modes at the
  + [`BuildProjectionOperatorsForIGWModes`](/classes/transforms/wvtransformboussinesq/buildprojectionoperatorsforigwmodes.html) Now go compute the appropriate number of modes at the
  + [`BuildProjectionOperatorsWithFreeSurface`](/classes/transforms/wvtransformboussinesq/buildprojectionoperatorswithfreesurface.html) Make these matrices invertible by adding the barotropic mode
  + [`BuildProjectionOperatorsWithRigidLid`](/classes/transforms/wvtransformboussinesq/buildprojectionoperatorswithrigidlid.html) 
  + [`K2unique`](/classes/transforms/wvtransformboussinesq/k2unique.html) unique squared-wavenumbers
  + [`K2uniqueK2Map`](/classes/transforms/wvtransformboussinesq/k2uniquek2map.html) cell array Nk in length. Each cell contains indices back to K2
  + [`N2`](/classes/transforms/wvtransformboussinesq/n2.html) 
  + [`N2Function`](/classes/transforms/wvtransformboussinesq/n2function.html) 
  + [`P0`](/classes/transforms/wvtransformboussinesq/p0.html) Preconditioner for F, size(P)=[1 1 Nj x 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
  + [`PF0`](/classes/transforms/wvtransformboussinesq/pf0.html) size(PF,PG)=[Nj x Nz x 1]
  + [`PF0inv`](/classes/transforms/wvtransformboussinesq/pf0inv.html) Geostrophic transformation matrices
  + [`PFinvInterp`](/classes/transforms/wvtransformboussinesq/pfinvinterp.html) 
  + [`PFpm`](/classes/transforms/wvtransformboussinesq/pfpm.html) size(PF,PG)=[Nj x Nz x Nk]
  + [`PFpmInv`](/classes/transforms/wvtransformboussinesq/pfpminv.html) IGW transformation matrices
  + [`Ppm`](/classes/transforms/wvtransformboussinesq/ppm.html) Preconditioner for F, size(P)=[Nj x Nk]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
  + [`Q0`](/classes/transforms/wvtransformboussinesq/q0.html) Preconditioner for G, size(Q)=[1 1 Nj x 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
  + [`QG0`](/classes/transforms/wvtransformboussinesq/qg0.html) 
  + [`QG0inv`](/classes/transforms/wvtransformboussinesq/qg0inv.html) 
  + [`QGinvInterp`](/classes/transforms/wvtransformboussinesq/qginvinterp.html) 
  + [`QGpm`](/classes/transforms/wvtransformboussinesq/qgpm.html) 
  + [`QGpmInv`](/classes/transforms/wvtransformboussinesq/qgpminv.html) 
  + [`QGwg`](/classes/transforms/wvtransformboussinesq/qgwg.html) size(PF,PG)=[Nj x Nj x Nk]
  + [`Qpm`](/classes/transforms/wvtransformboussinesq/qpm.html) Preconditioner for G, size(Q)=[Nj x Nk]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
  + [`WVTransformBoussinesq`](/classes/transforms/wvtransformboussinesq/wvtransformboussinesq.html) if all of these things are set initially (presumably read
  + [`Wklz`](/classes/transforms/wvtransformboussinesq/wklz.html) 
  + [`Wzkl`](/classes/transforms/wvtransformboussinesq/wzkl.html) 
  + [`buildInterpolationProjectionOperators`](/classes/transforms/wvtransformboussinesq/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes/transforms/wvtransformboussinesq/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`dLnN2`](/classes/transforms/wvtransformboussinesq/dlnn2.html) 
  + [`dLnN2Function`](/classes/transforms/wvtransformboussinesq/dlnn2function.html) 
  + [`dftBuffer`](/classes/transforms/wvtransformboussinesq/dftbuffer.html) 
  + [`dftConjugateIndex`](/classes/transforms/wvtransformboussinesq/dftconjugateindex.html) 
  + [`dftPrimaryIndex`](/classes/transforms/wvtransformboussinesq/dftprimaryindex.html) 
  + [`diffZF`](/classes/transforms/wvtransformboussinesq/diffzf.html) 
  + [`diffZG`](/classes/transforms/wvtransformboussinesq/diffzg.html) 
  + [`h_0`](/classes/transforms/wvtransformboussinesq/h_0.html) [Nj 1]
  + [`h_pm`](/classes/transforms/wvtransformboussinesq/h_pm.html) [Nj Nkl]
  + [`iK2unique`](/classes/transforms/wvtransformboussinesq/ik2unique.html) map from 2-dim K2, to 1-dim K2unique
  + [`iOmega`](/classes/transforms/wvtransformboussinesq/iomega.html) 
  + [`internalModes`](/classes/transforms/wvtransformboussinesq/internalmodes.html) 
  + [`isHydrostatic`](/classes/transforms/wvtransformboussinesq/ishydrostatic.html) 
  + [`isequal`](/classes/transforms/wvtransformboussinesq/isequal.html) 
  + [`nK2unique`](/classes/transforms/wvtransformboussinesq/nk2unique.html) number of unique squared-wavenumbers
  + [`rhoFunction`](/classes/transforms/wvtransformboussinesq/rhofunction.html) function handles
  + [`rhobar`](/classes/transforms/wvtransformboussinesq/rhobar.html) on the z-grid, size(N2) = [length(z) 1];
  + [`transformFromSpatialDomainWithFg`](/classes/transforms/wvtransformboussinesq/transformfromspatialdomainwithfg.html) 
  + [`transformFromSpatialDomainWithFio`](/classes/transforms/wvtransformboussinesq/transformfromspatialdomainwithfio.html) Required for transformUVEtaToWaveVortex
  + [`transformFromSpatialDomainWithFourier`](/classes/transforms/wvtransformboussinesq/transformfromspatialdomainwithfourier.html) self.dftBuffer = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
  + [`transformFromSpatialDomainWithGg`](/classes/transforms/wvtransformboussinesq/transformfromspatialdomainwithgg.html) 
  + [`transformToSpatialDomainWithF`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithf.html) Required for transformWaveVortexToUVEta
  + [`transformToSpatialDomainWithFInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithFg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfg.html) arguments
  + [`transformToSpatialDomainWithFw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfw.html) 
  + [`transformToSpatialDomainWithG`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithg.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithginterp.html) 
  + [`transformToSpatialDomainWithGg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgg.html) arguments
  + [`transformToSpatialDomainWithGw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgw.html) 
  + [`transformWithG_wg`](/classes/transforms/wvtransformboussinesq/transformwithg_wg.html) 
  + [`uMaxA0`](/classes/transforms/wvtransformboussinesq/umaxa0.html) 
  + [`uMaxGNormRatioForWave`](/classes/transforms/wvtransformboussinesq/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`wvBuffer`](/classes/transforms/wvtransformboussinesq/wvbuffer.html) 
  + [`wvConjugateIndex`](/classes/transforms/wvtransformboussinesq/wvconjugateindex.html) 
  + [`zInterp`](/classes/transforms/wvtransformboussinesq/zinterp.html) 
+ Initialization (Static)
  + [`waveVortexTransformFromFile`](/classes/transforms/wvtransformboussinesq/wavevortextransformfromfile.html) Initialize a WVTransformHydrostatic instance from an existing file
+ Initialization
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformboussinesq/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Write to file
  + [`writeToFile`](/classes/transforms/wvtransformboussinesq/writetofile.html) Output the `WVTransform` to file.


---