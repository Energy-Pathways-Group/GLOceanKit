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
+ Stratification
  + [`N2`](/classes/transforms/wvtransformboussinesq/n2.html) buoyancy frequency of the no-motion density
  + [`dLnN2`](/classes/transforms/wvtransformboussinesq/dlnn2.html) d/dz ln N2
  + [`rho_nm`](/classes/transforms/wvtransformboussinesq/rho_nm.html) no-motion density
  + [`verticalModes`](/classes/transforms/wvtransformboussinesq/verticalmodes.html) instance of the InternalModes class
  + Vertical modes
    + [`FMatrix`](/classes/transforms/wvtransformboussinesq/fmatrix.html) transformation matrix $$F_g$$
    + [`FinvMatrix`](/classes/transforms/wvtransformboussinesq/finvmatrix.html) transformation matrix $$F_g^{-1}$$
    + [`GMatrix`](/classes/transforms/wvtransformboussinesq/gmatrix.html) transformation matrix $$G_g$$
    + [`GinvMatrix`](/classes/transforms/wvtransformboussinesq/ginvmatrix.html) transformation matrix $$G_g^{-1}$$
+ Initial conditions
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/transforms/wvtransformboussinesq/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/transforms/wvtransformboussinesq/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/transforms/wvtransformboussinesq/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/transforms/wvtransformboussinesq/setinertialmotions.html) set inertial motions
+ Operations
  + Differentiation
    + [`diffZF`](/classes/transforms/wvtransformboussinesq/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes/transforms/wvtransformboussinesq/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
+ Other
  + [`WVTransformBoussinesq`](/classes/transforms/wvtransformboussinesq/wvtransformboussinesq.html) First we need to initialize the WVStratifiedFlow.
  + [`buildInterpolationProjectionOperators`](/classes/transforms/wvtransformboussinesq/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes/transforms/wvtransformboussinesq/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`buildVerticalModeProjectionOperators`](/classes/transforms/wvtransformboussinesq/buildverticalmodeprojectionoperators.html) 
  + [`nK2unique`](/classes/transforms/wvtransformboussinesq/nk2unique.html) number of unique squared-wavenumbers
  + [`transformToSpatialDomainWithFInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithFg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfg.html) arguments
  + [`transformToSpatialDomainWithFw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfw.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithginterp.html) 
  + [`transformToSpatialDomainWithGg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgg.html) arguments
  + [`transformToSpatialDomainWithGw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgw.html) 
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformboussinesq/wavevortextransformwithresolution.html) 


---