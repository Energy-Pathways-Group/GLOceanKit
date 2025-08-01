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
+ Primary flow components
  + [`geostrophicComponent`](/classes/transforms/wvtransformhydrostatic/geostrophiccomponent.html) returns the geostrophic flow component
  + [`waveComponent`](/classes/transforms/wvtransformhydrostatic/wavecomponent.html) returns the internal gravity wave flow component
  + [`inertialComponent`](/classes/transforms/wvtransformhydrostatic/inertialcomponent.html) returns the inertial oscillation flow component
  + [`mdaComponent`](/classes/transforms/wvtransformhydrostatic/mdacomponent.html) returns the mean density anomaly component
+ Stratification
  + [`rho_nm`](/classes/transforms/wvtransformhydrostatic/rho_nm.html) $$\rho_\textrm{nm}(z)$$, no-motion density
  + [`N2`](/classes/transforms/wvtransformhydrostatic/n2.html) $$N^2(z)$$, squared buoyancy frequency of the no-motion density, $$N^2\equiv - \frac{g}{\rho_0} \frac{\partial \rho_\textrm{nm}}{\partial z}$$
  + [`dLnN2`](/classes/transforms/wvtransformhydrostatic/dlnn2.html) $$\frac{\partial \ln N^2}{\partial z}$$, vertical variation of the log of the squared buoyancy frequency
  + [`verticalModes`](/classes/transforms/wvtransformhydrostatic/verticalmodes.html) instance of the InternalModes class
  + [`effectiveVerticalGridResolution`](/classes/transforms/wvtransformhydrostatic/effectiveverticalgridresolution.html) returns the effective vertical grid resolution in meters
  + Vertical modes
    + [`FMatrix`](/classes/transforms/wvtransformhydrostatic/fmatrix.html) transformation matrix $$F_g$$
    + [`FinvMatrix`](/classes/transforms/wvtransformhydrostatic/finvmatrix.html) transformation matrix $$F_g^{-1}$$
    + [`GMatrix`](/classes/transforms/wvtransformhydrostatic/gmatrix.html) transformation matrix $$G_g$$
    + [`GinvMatrix`](/classes/transforms/wvtransformhydrostatic/ginvmatrix.html) transformation matrix $$G_g^{-1}$$
  + Validation
    + [`isDensityInValidRange`](/classes/transforms/wvtransformhydrostatic/isdensityinvalidrange.html) checks if the density field is a valid adiabatic re-arrangement of the base state
+ Initial conditions
  + Waves
    + [`addGMSpectrum`](/classes/transforms/wvtransformhydrostatic/addgmspectrum.html) add waves following a Garrett-Munk spectrum
    + [`addWaveModes`](/classes/transforms/wvtransformhydrostatic/addwavemodes.html) add amplitudes of the given wave modes
    + [`addWavesWithFrequencySpectrum`](/classes/transforms/wvtransformhydrostatic/addwaveswithfrequencyspectrum.html) add waves with a specified frequency spectrum
    + [`initWavesWithFrequencySpectrum`](/classes/transforms/wvtransformhydrostatic/initwaveswithfrequencyspectrum.html) initialize with waves of a specified frequency spectrum
    + [`initWithAlternativeSpectrum`](/classes/transforms/wvtransformhydrostatic/initwithalternativespectrum.html) initialize with an alternative formulation of the GM spectrum in the wavenumber domain.
    + [`initWithGMSpectrum`](/classes/transforms/wvtransformhydrostatic/initwithgmspectrum.html) initialize the wave field following a Garrett-Munk spectrum
    + [`initWithWaveModes`](/classes/transforms/wvtransformhydrostatic/initwithwavemodes.html) initialize with the given wave modes
    + [`removeAllWaves`](/classes/transforms/wvtransformhydrostatic/removeallwaves.html) removes all wave from the model, including inertial oscillations
    + [`setWaveModes`](/classes/transforms/wvtransformhydrostatic/setwavemodes.html) set amplitudes of the given wave modes
  + Geostrophic Motions
    + [`initWithGeostrophicStreamfunction`](/classes/transforms/wvtransformhydrostatic/initwithgeostrophicstreamfunction.html) initialize with a geostrophic streamfunction
    + [`setGeostrophicStreamfunction`](/classes/transforms/wvtransformhydrostatic/setgeostrophicstreamfunction.html) set a geostrophic streamfunction
    + [`addGeostrophicStreamfunction`](/classes/transforms/wvtransformhydrostatic/addgeostrophicstreamfunction.html) add a geostrophic streamfunction to existing geostrophic motions
    + [`setGeostrophicModes`](/classes/transforms/wvtransformhydrostatic/setgeostrophicmodes.html) set amplitudes of the given geostrophic modes
    + [`addGeostrophicModes`](/classes/transforms/wvtransformhydrostatic/addgeostrophicmodes.html) add amplitudes of the given geostrophic modes
    + [`removeAllGeostrophicMotions`](/classes/transforms/wvtransformhydrostatic/removeallgeostrophicmotions.html) remove all geostrophic motions
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/transforms/wvtransformhydrostatic/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/transforms/wvtransformhydrostatic/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/transforms/wvtransformhydrostatic/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/transforms/wvtransformhydrostatic/setinertialmotions.html) set inertial motions
  + Mean density anomaly
    + [`addMeanDensityAnomaly`](/classes/transforms/wvtransformhydrostatic/addmeandensityanomaly.html) add inertial motions to existing inertial motions
    + [`initWithMeanDensityAnomaly`](/classes/transforms/wvtransformhydrostatic/initwithmeandensityanomaly.html) initialize with inertial motions
    + [`removeAllMeanDensityAnomaly`](/classes/transforms/wvtransformhydrostatic/removeallmeandensityanomaly.html) remove all mean density anomalies
    + [`setMeanDensityAnomaly`](/classes/transforms/wvtransformhydrostatic/setmeandensityanomaly.html) set inertial motions
+ Energetics of flow components
  + [`geostrophicEnergy`](/classes/transforms/wvtransformhydrostatic/geostrophicenergy.html) total energy, geostrophic
  + [`waveEnergy`](/classes/transforms/wvtransformhydrostatic/waveenergy.html) total energy, waves
  + [`inertialEnergy`](/classes/transforms/wvtransformhydrostatic/inertialenergy.html) total energy, inertial oscillations
  + [`mdaEnergy`](/classes/transforms/wvtransformhydrostatic/mdaenergy.html) total energy, mean density anomaly
+ Operations
  + Differentiation
    + [`diffZF`](/classes/transforms/wvtransformhydrostatic/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes/transforms/wvtransformhydrostatic/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
+ Energetics
  + [`geostrophicKineticEnergy`](/classes/transforms/wvtransformhydrostatic/geostrophickineticenergy.html) kinetic energy of the geostrophic flow
  + [`geostrophicPotentialEnergy`](/classes/transforms/wvtransformhydrostatic/geostrophicpotentialenergy.html) potential energy of the geostrophic flow
+ Other
  + [`N2Function`](/classes/transforms/wvtransformhydrostatic/n2function.html) 
  + [`dLnN2Function`](/classes/transforms/wvtransformhydrostatic/dlnn2function.html) 
  + [`maxFg`](/classes/transforms/wvtransformhydrostatic/maxfg.html) 
  + [`maxFw`](/classes/transforms/wvtransformhydrostatic/maxfw.html) 
  + [`rhoFunction`](/classes/transforms/wvtransformhydrostatic/rhofunction.html) function handles
  + [`throwErrorIfDensityViolation`](/classes/transforms/wvtransformhydrostatic/throwerrorifdensityviolation.html) checks if the proposed coefficients are a valid adiabatic re-arrangement of the base state


---