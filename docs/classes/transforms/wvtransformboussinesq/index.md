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
  + [`rho_nm`](/classes/transforms/wvtransformboussinesq/rho_nm.html) $$\rho_\textrm{nm}(z)$$, no-motion density
  + [`N2`](/classes/transforms/wvtransformboussinesq/n2.html) $$N^2(z)$$, squared buoyancy frequency of the no-motion density, $$N^2\equiv - \frac{g}{\rho_0} \frac{\partial \rho_\textrm{nm}}{\partial z}$$
  + [`dLnN2`](/classes/transforms/wvtransformboussinesq/dlnn2.html) $$\frac{\partial \ln N^2}{\partial z}$$, vertical variation of the log of the squared buoyancy frequency
  + [`verticalModes`](/classes/transforms/wvtransformboussinesq/verticalmodes.html) instance of the InternalModes class
  + [`effectiveVerticalGridResolution`](/classes/transforms/wvtransformboussinesq/effectiveverticalgridresolution.html) returns the effective vertical grid resolution in meters
  + Vertical modes
    + [`FMatrix`](/classes/transforms/wvtransformboussinesq/fmatrix.html) transformation matrix $$F_g$$
    + [`FinvMatrix`](/classes/transforms/wvtransformboussinesq/finvmatrix.html) transformation matrix $$F_g^{-1}$$
    + [`GMatrix`](/classes/transforms/wvtransformboussinesq/gmatrix.html) transformation matrix $$G_g$$
    + [`GinvMatrix`](/classes/transforms/wvtransformboussinesq/ginvmatrix.html) transformation matrix $$G_g^{-1}$$
  + Validation
    + [`isDensityInValidRange`](/classes/transforms/wvtransformboussinesq/isdensityinvalidrange.html) checks if the density field is a valid adiabatic re-arrangement of the base state
+ Initial conditions
  + Geostrophic Motions
    + [`initWithGeostrophicStreamfunction`](/classes/transforms/wvtransformboussinesq/initwithgeostrophicstreamfunction.html) initialize with a geostrophic streamfunction
    + [`setGeostrophicStreamfunction`](/classes/transforms/wvtransformboussinesq/setgeostrophicstreamfunction.html) set a geostrophic streamfunction
    + [`addGeostrophicStreamfunction`](/classes/transforms/wvtransformboussinesq/addgeostrophicstreamfunction.html) add a geostrophic streamfunction to existing geostrophic motions
    + [`setGeostrophicModes`](/classes/transforms/wvtransformboussinesq/setgeostrophicmodes.html) set amplitudes of the given geostrophic modes
    + [`addGeostrophicModes`](/classes/transforms/wvtransformboussinesq/addgeostrophicmodes.html) add amplitudes of the given geostrophic modes
    + [`removeAllGeostrophicMotions`](/classes/transforms/wvtransformboussinesq/removeallgeostrophicmotions.html) remove all geostrophic motions
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/transforms/wvtransformboussinesq/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/transforms/wvtransformboussinesq/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/transforms/wvtransformboussinesq/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/transforms/wvtransformboussinesq/setinertialmotions.html) set inertial motions
  + Mean density anomaly
    + [`addMeanDensityAnomaly`](/classes/transforms/wvtransformboussinesq/addmeandensityanomaly.html) add inertial motions to existing inertial motions
    + [`initWithMeanDensityAnomaly`](/classes/transforms/wvtransformboussinesq/initwithmeandensityanomaly.html) initialize with inertial motions
    + [`removeAllMeanDensityAnomaly`](/classes/transforms/wvtransformboussinesq/removeallmeandensityanomaly.html) remove all mean density anomalies
    + [`setMeanDensityAnomaly`](/classes/transforms/wvtransformboussinesq/setmeandensityanomaly.html) set inertial motions
  + Waves
    + [`addWaveModes`](/classes/transforms/wvtransformboussinesq/addwavemodes.html) add amplitudes of the given wave modes
    + [`initWithHorizontalWaveNumberSpectrum`](/classes/transforms/wvtransformboussinesq/initwithhorizontalwavenumberspectrum.html) initialize with a Alternative Interal Wave Spectrum in
    + [`initWithWaveModes`](/classes/transforms/wvtransformboussinesq/initwithwavemodes.html) initialize with the given wave modes
    + [`removeAllWaves`](/classes/transforms/wvtransformboussinesq/removeallwaves.html) removes all wave from the model, including inertial oscillations
    + [`setWaveModes`](/classes/transforms/wvtransformboussinesq/setwavemodes.html) set amplitudes of the given wave modes
+ Operations
  + Differentiation
    + [`diffZF`](/classes/transforms/wvtransformboussinesq/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes/transforms/wvtransformboussinesq/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
+ Primary flow components
  + [`geostrophicComponent`](/classes/transforms/wvtransformboussinesq/geostrophiccomponent.html) returns the geostrophic flow component
  + [`waveComponent`](/classes/transforms/wvtransformboussinesq/wavecomponent.html) returns the internal gravity wave flow component
  + [`inertialComponent`](/classes/transforms/wvtransformboussinesq/inertialcomponent.html) returns the inertial oscillation flow component
  + [`mdaComponent`](/classes/transforms/wvtransformboussinesq/mdacomponent.html) returns the mean density anomaly component
+ Energetics of flow components
  + [`geostrophicEnergy`](/classes/transforms/wvtransformboussinesq/geostrophicenergy.html) total energy, geostrophic
  + [`waveEnergy`](/classes/transforms/wvtransformboussinesq/waveenergy.html) total energy, waves
  + [`inertialEnergy`](/classes/transforms/wvtransformboussinesq/inertialenergy.html) total energy, inertial oscillations
  + [`mdaEnergy`](/classes/transforms/wvtransformboussinesq/mdaenergy.html) total energy, mean density anomaly
+ Energetics
  + [`geostrophicKineticEnergy`](/classes/transforms/wvtransformboussinesq/geostrophickineticenergy.html) kinetic energy of the geostrophic flow
  + [`geostrophicPotentialEnergy`](/classes/transforms/wvtransformboussinesq/geostrophicpotentialenergy.html) potential energy of the geostrophic flow
+ Other
  + [`N2Function`](/classes/transforms/wvtransformboussinesq/n2function.html) 
  + [`WVTransformBoussinesq`](/classes/transforms/wvtransformboussinesq/wvtransformboussinesq.html) First we need to initialize the WVStratifiedFlow.
  + [`buildInterpolationProjectionOperators`](/classes/transforms/wvtransformboussinesq/buildinterpolationprojectionoperators.html) 
  + [`buildInterpolationProjectionOperatorsForGrid`](/classes/transforms/wvtransformboussinesq/buildinterpolationprojectionoperatorsforgrid.html) 
  + [`buildVerticalModeProjectionOperators`](/classes/transforms/wvtransformboussinesq/buildverticalmodeprojectionoperators.html) 
  + [`dLnN2Function`](/classes/transforms/wvtransformboussinesq/dlnn2function.html) 
  + [`initWithAlternativeSpectrum`](/classes/transforms/wvtransformboussinesq/initwithalternativespectrum.html) 
  + [`initWithFrequencySpectrum`](/classes/transforms/wvtransformboussinesq/initwithfrequencyspectrum.html) 
  + [`initWithGMSpectrum`](/classes/transforms/wvtransformboussinesq/initwithgmspectrum.html) 
  + [`maxFg`](/classes/transforms/wvtransformboussinesq/maxfg.html) 
  + [`maxFw`](/classes/transforms/wvtransformboussinesq/maxfw.html) 
  + [`nK2unique`](/classes/transforms/wvtransformboussinesq/nk2unique.html) number of unique squared-wavenumbers
  + [`rhoFunction`](/classes/transforms/wvtransformboussinesq/rhofunction.html) function handles
  + [`throwErrorIfDensityViolation`](/classes/transforms/wvtransformboussinesq/throwerrorifdensityviolation.html) checks if the proposed coefficients are a valid adiabatic re-arrangement of the base state
  + [`transformToSpatialDomainWithFInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfinterp.html) 
  + [`transformToSpatialDomainWithFg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfg.html) arguments
  + [`transformToSpatialDomainWithFw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithfw.html) 
  + [`transformToSpatialDomainWithGInterp`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithginterp.html) 
  + [`transformToSpatialDomainWithGg`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgg.html) arguments
  + [`transformToSpatialDomainWithGw`](/classes/transforms/wvtransformboussinesq/transformtospatialdomainwithgw.html) 


---