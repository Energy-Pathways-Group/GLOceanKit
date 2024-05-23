---
layout: default
title: WVTransformHydrostatic
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
nav_order: 2
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
+ Stratification
  + [`N2`](/classes/transforms/wvtransformhydrostatic/n2.html) buoyancy frequency of the no-motion density
  + [`dLnN2`](/classes/transforms/wvtransformhydrostatic/dlnn2.html) d/dz ln N2
  + [`rho_nm`](/classes/transforms/wvtransformhydrostatic/rho_nm.html) no-motion density
  + [`verticalModes`](/classes/transforms/wvtransformhydrostatic/verticalmodes.html) instance of the InternalModes class
  + Vertical modes
    + [`FMatrix`](/classes/transforms/wvtransformhydrostatic/fmatrix.html) transformation matrix $$F_g$$
    + [`FinvMatrix`](/classes/transforms/wvtransformhydrostatic/finvmatrix.html) transformation matrix $$F_g^{-1}$$
    + [`GMatrix`](/classes/transforms/wvtransformhydrostatic/gmatrix.html) transformation matrix $$G_g$$
    + [`GinvMatrix`](/classes/transforms/wvtransformhydrostatic/ginvmatrix.html) transformation matrix $$G_g^{-1}$$
+ Initial conditions
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/transforms/wvtransformhydrostatic/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/transforms/wvtransformhydrostatic/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/transforms/wvtransformhydrostatic/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/transforms/wvtransformhydrostatic/setinertialmotions.html) set inertial motions
+ Operations
  + Differentiation
    + [`diffZF`](/classes/transforms/wvtransformhydrostatic/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes/transforms/wvtransformhydrostatic/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
+ Other
  + [`waveVortexTransformWithResolution`](/classes/transforms/wvtransformhydrostatic/wavevortextransformwithresolution.html) 


---