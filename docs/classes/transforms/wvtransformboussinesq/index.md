---
layout: default
title: WVTransformBoussinesq
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
nav_order: 1
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
+ Initial conditions
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/transforms/wvtransformboussinesq/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/transforms/wvtransformboussinesq/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/transforms/wvtransformboussinesq/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/transforms/wvtransformboussinesq/setinertialmotions.html) set inertial motions
+ Other
  + [`Aklz`](/classes/transforms/wvtransformboussinesq/aklz.html) 
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
  + [`buildVerticalModeProjectionOperators`](/classes/transforms/wvtransformboussinesq/buildverticalmodeprojectionoperators.html) 
  + [`dLnN2`](/classes/transforms/wvtransformboussinesq/dlnn2.html) 
  + [`dLnN2Function`](/classes/transforms/wvtransformboussinesq/dlnn2function.html) 
  + [`dftBuffer`](/classes/transforms/wvtransformboussinesq/dftbuffer.html) 
  + [`dftConjugateIndex`](/classes/transforms/wvtransformboussinesq/dftconjugateindex.html) 
  + [`dftPrimaryIndex`](/classes/transforms/wvtransformboussinesq/dftprimaryindex.html) 
  + [`iK2unique`](/classes/transforms/wvtransformboussinesq/ik2unique.html) map from 2-dim K2, to 1-dim K2unique
  + [`nK2unique`](/classes/transforms/wvtransformboussinesq/nk2unique.html) number of unique squared-wavenumbers
  + [`rhoFunction`](/classes/transforms/wvtransformboussinesq/rhofunction.html) function handles
  + [`rhobar`](/classes/transforms/wvtransformboussinesq/rhobar.html) on the z-grid, size(N2) = [length(z) 1];
  + [`transformToSpatialDomainWithFInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithFg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfg.html) arguments
  + [`transformToSpatialDomainWithFw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfw.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithginterp.html) 
  + [`transformToSpatialDomainWithGg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgg.html) arguments
  + [`transformToSpatialDomainWithGw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgw.html) 
  + [`verticalModes`](/classes/transforms/wvtransformboussinesq/verticalmodes.html) 
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformboussinesq/wavevortextransformwithresolution.html) 
  + [`wvBuffer`](/classes/transforms/wvtransformboussinesq/wvbuffer.html) 
  + [`wvConjugateIndex`](/classes/transforms/wvtransformboussinesq/wvconjugateindex.html) 
  + [`zInterp`](/classes/transforms/wvtransformboussinesq/zinterp.html) 


---