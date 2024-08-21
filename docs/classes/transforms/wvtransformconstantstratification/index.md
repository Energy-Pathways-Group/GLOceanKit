---
layout: default
title: WVTransformConstantStratification
has_children: false
has_toc: false
mathjax: true
parent: Transforms
grand_parent: Class documentation
nav_order: 3
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
  wvt = WVTransformConstantStratification([100e3, 100e3, 1300],[64, 64, 65], N0=N0,latitude=30);
  ```
 
   
  


## Topics
+ Initialization
  + [`WVTransformConstantStratification`](/classes/transforms/wvtransformconstantstratification/wvtransformconstantstratification.html) initialze a wave-vortex transform with constant stratification
+ Stratification
  + [`rho_nm`](/classes/transforms/wvtransformconstantstratification/rho_nm.html) $$\rho_\textrm{nm}(z)$$, no-motion density
  + [`N2`](/classes/transforms/wvtransformconstantstratification/n2.html) $$N^2(z)$$, squared buoyancy frequency of the no-motion density, $$N^2\equiv - \frac{g}{\rho_0} \frac{\partial \rho_\textrm{nm}}{\partial z}$$
  + [`dLnN2`](/classes/transforms/wvtransformconstantstratification/dlnn2.html) $$\frac{\partial \ln N^2}{\partial z}$$, vertical variation of the log of the squared buoyancy frequency
  + [`verticalModes`](/classes/transforms/wvtransformconstantstratification/verticalmodes.html) instance of the InternalModes class
  + [`effectiveVerticalGridResolution`](/classes/transforms/wvtransformconstantstratification/effectiveverticalgridresolution.html) returns the effective vertical grid resolution in meters
  + Vertical modes
    + [`FMatrix`](/classes/transforms/wvtransformconstantstratification/fmatrix.html) transformation matrix $$F_g$$
    + [`FinvMatrix`](/classes/transforms/wvtransformconstantstratification/finvmatrix.html) transformation matrix $$F_g^{-1}$$
    + [`GMatrix`](/classes/transforms/wvtransformconstantstratification/gmatrix.html) transformation matrix $$G_g$$
    + [`GinvMatrix`](/classes/transforms/wvtransformconstantstratification/ginvmatrix.html) transformation matrix $$G_g^{-1}$$
  + Validation
    + [`isDensityInValidRange`](/classes/transforms/wvtransformconstantstratification/isdensityinvalidrange.html) checks if the density field is a valid adiabatic re-arrangement of the base state
+ Initial conditions
  + Geostrophic Motions
    + [`initWithGeostrophicStreamfunction`](/classes/transforms/wvtransformconstantstratification/initwithgeostrophicstreamfunction.html) initialize with a geostrophic streamfunction
    + [`setGeostrophicStreamfunction`](/classes/transforms/wvtransformconstantstratification/setgeostrophicstreamfunction.html) set a geostrophic streamfunction
    + [`addGeostrophicStreamfunction`](/classes/transforms/wvtransformconstantstratification/addgeostrophicstreamfunction.html) add a geostrophic streamfunction to existing geostrophic motions
    + [`setGeostrophicModes`](/classes/transforms/wvtransformconstantstratification/setgeostrophicmodes.html) set amplitudes of the given geostrophic modes
    + [`addGeostrophicModes`](/classes/transforms/wvtransformconstantstratification/addgeostrophicmodes.html) add amplitudes of the given geostrophic modes
    + [`removeAllGeostrophicMotions`](/classes/transforms/wvtransformconstantstratification/removeallgeostrophicmotions.html) remove all geostrophic motions
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/transforms/wvtransformconstantstratification/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/transforms/wvtransformconstantstratification/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/transforms/wvtransformconstantstratification/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/transforms/wvtransformconstantstratification/setinertialmotions.html) set inertial motions
  + Mean density anomaly
    + [`addMeanDensityAnomaly`](/classes/transforms/wvtransformconstantstratification/addmeandensityanomaly.html) add inertial motions to existing inertial motions
    + [`initWithMeanDensityAnomaly`](/classes/transforms/wvtransformconstantstratification/initwithmeandensityanomaly.html) initialize with inertial motions
    + [`removeAllMeanDensityAnomaly`](/classes/transforms/wvtransformconstantstratification/removeallmeandensityanomaly.html) remove all mean density anomalies
    + [`setMeanDensityAnomaly`](/classes/transforms/wvtransformconstantstratification/setmeandensityanomaly.html) set inertial motions
  + Waves
    + [`addWaveModes`](/classes/transforms/wvtransformconstantstratification/addwavemodes.html) add amplitudes of the given wave modes
    + [`initWithHorizontalWaveNumberSpectrum`](/classes/transforms/wvtransformconstantstratification/initwithhorizontalwavenumberspectrum.html) initialize with a Alternative Interal Wave Spectrum in
    + [`initWithWaveModes`](/classes/transforms/wvtransformconstantstratification/initwithwavemodes.html) initialize with the given wave modes
    + [`removeAllWaves`](/classes/transforms/wvtransformconstantstratification/removeallwaves.html) removes all wave from the model, including inertial oscillations
    + [`setWaveModes`](/classes/transforms/wvtransformconstantstratification/setwavemodes.html) set amplitudes of the given wave modes
+ Operations
  + Differentiation
    + [`diffZF`](/classes/transforms/wvtransformconstantstratification/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes/transforms/wvtransformconstantstratification/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
+ Primary flow components
  + [`geostrophicComponent`](/classes/transforms/wvtransformconstantstratification/geostrophiccomponent.html) returns the geostrophic flow component
  + [`waveComponent`](/classes/transforms/wvtransformconstantstratification/wavecomponent.html) returns the internal gravity wave flow component
  + [`inertialComponent`](/classes/transforms/wvtransformconstantstratification/inertialcomponent.html) returns the inertial oscillation flow component
  + [`mdaComponent`](/classes/transforms/wvtransformconstantstratification/mdacomponent.html) returns the mean density anomaly component
+ Energetics of flow components
  + [`geostrophicEnergy`](/classes/transforms/wvtransformconstantstratification/geostrophicenergy.html) total energy, geostrophic
  + [`waveEnergy`](/classes/transforms/wvtransformconstantstratification/waveenergy.html) total energy, waves
  + [`inertialEnergy`](/classes/transforms/wvtransformconstantstratification/inertialenergy.html) total energy, inertial oscillations
  + [`mdaEnergy`](/classes/transforms/wvtransformconstantstratification/mdaenergy.html) total energy, mean density anomaly
+ Energetics
  + [`geostrophicKineticEnergy`](/classes/transforms/wvtransformconstantstratification/geostrophickineticenergy.html) kinetic energy of the geostrophic flow
  + [`geostrophicPotentialEnergy`](/classes/transforms/wvtransformconstantstratification/geostrophicpotentialenergy.html) potential energy of the geostrophic flow
+ Other
  + [`CosineTransformBackMatrix`](/classes/transforms/wvtransformconstantstratification/cosinetransformbackmatrix.html) Discrete Cosine Transform (DCT-I) matrix
  + [`CosineTransformForwardMatrix`](/classes/transforms/wvtransformconstantstratification/cosinetransformforwardmatrix.html) Discrete Cosine Transform (DCT-I) matrix
  + [`N0`](/classes/transforms/wvtransformconstantstratification/n0.html) 
  + [`N2AtDepth`](/classes/transforms/wvtransformconstantstratification/n2atdepth.html) 
  + [`N2Function`](/classes/transforms/wvtransformconstantstratification/n2function.html) 
  + [`PlaceParticlesOnIsopycnal`](/classes/transforms/wvtransformconstantstratification/placeparticlesonisopycnal.html) MAS 1/10/18 - added intext ('int' or 'both') to give option of using int vs. int+ext fields for rho_prime
  + [`RhoBarAtDepth`](/classes/transforms/wvtransformconstantstratification/rhobaratdepth.html) 
  + [`SineTransformBackMatrix`](/classes/transforms/wvtransformconstantstratification/sinetransformbackmatrix.html) CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix
  + [`SineTransformForwardMatrix`](/classes/transforms/wvtransformconstantstratification/sinetransformforwardmatrix.html) CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix
  + [`buildVerticalModeProjectionOperators`](/classes/transforms/wvtransformconstantstratification/buildverticalmodeprojectionoperators.html) We renormalization the transformation matrices to directly
  + [`dLnN2Function`](/classes/transforms/wvtransformconstantstratification/dlnn2function.html) 
  + [`initWithAlternativeSpectrum`](/classes/transforms/wvtransformconstantstratification/initwithalternativespectrum.html) 
  + [`initWithFrequencySpectrum`](/classes/transforms/wvtransformconstantstratification/initwithfrequencyspectrum.html) 
  + [`initWithGMSpectrum`](/classes/transforms/wvtransformconstantstratification/initwithgmspectrum.html) 
  + [`maxFg`](/classes/transforms/wvtransformconstantstratification/maxfg.html) 
  + [`maxFw`](/classes/transforms/wvtransformconstantstratification/maxfw.html) 
  + [`rhoFunction`](/classes/transforms/wvtransformconstantstratification/rhofunction.html) function handles
  + [`throwErrorIfDensityViolation`](/classes/transforms/wvtransformconstantstratification/throwerrorifdensityviolation.html) checks if the proposed coefficients are a valid adiabatic re-arrangement of the base state
  + [`transformFromSpatialDomainWithF`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithf.html) 
  + [`transformFromSpatialDomainWithG`](/classes/transforms/wvtransformconstantstratification/transformfromspatialdomainwithg.html) 


---