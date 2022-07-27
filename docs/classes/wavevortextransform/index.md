---
layout: default
title: WaveVortexTransform
parent: Classes
has_children: false
has_toc: false
mathjax: true
---

#  WaveVortexTransform

Represents the state of the ocean in terms of energetically orthogonal wave and geostrophic (vortex) solutions

## Overview
 
 
  The WaveVortexTransform subclasses encapsulate data representing the
  state of the ocean at a given instant in time. What makes the
  WaveVortexTransform subclasses special is that the state of the ocean
  is represented as energetically independent waves and geostrophic
  motions (vortices). These classes can be queried for any ocean state
  variable including $$u$$, $$v$$, $$w$$, $$\rho$$, $$p$$, but also
  Ertel PV, relative vorticity, or custom defined state variables.
 
  The WaveVortexTransform is an abstract class and as such you must
  instatiate one of the concrete subclasses,
 
  + `WaveVortexTransformConstantStratification`
  + `WaveVortexTransformHydrostatic`
  + `WaveVortexTransformSingleMode`
 
                        


## Topics
+ Initialization
  + [`waveVortexTransformFromFile`](/classes/wavevortextransform/wavevortextransformfromfile.html) Initialize a WaveVortexTransform instance from an existing file
  + [`waveVortexTransformWithDoubleResolution`](/classes/wavevortextransform/wavevortextransformwithdoubleresolution.html) create a new WaveVortexTransform with double resolution
  + [`waveVortexTransformWithResolution`](/classes/wavevortextransform/wavevortextransformwithresolution.html) create a new WaveVortexTransform with increased resolution
+ Domain attributes
  + [`Omega`](/classes/wavevortextransform/omega.html) frequency of oscillation of the linear waves
  + [`f`](/classes/wavevortextransform/f.html) Coriolis parameter
  + [`g`](/classes/wavevortextransform/g.html) gravity of Earth
  + [`inertialPeriod`](/classes/wavevortextransform/inertialperiod.html) inertial period
  + [`isBarotropic`](/classes/wavevortextransform/isbarotropic.html) Boolean indicating whether there is a single (equivalent barotropic) mode
  + [`kRadial`](/classes/wavevortextransform/kradial.html) isotropic wavenumber dimension
  + [`latitude`](/classes/wavevortextransform/latitude.html) latitude of the simulation
  + [`t`](/classes/wavevortextransform/t.html) time of observations
  + [`t0`](/classes/wavevortextransform/t0.html) reference time of Ap, Am, A0
  + Grid
    + [`J`](/classes/wavevortextransform/j.html) j-coordinate matrix
    + [`K`](/classes/wavevortextransform/k.html) k-coordinate matrix
    + [`Kh`](/classes/wavevortextransform/kh.html) horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$
    + [`L`](/classes/wavevortextransform/l.html) l-coordinate matrix
    + [`Lx`](/classes/wavevortextransform/lx.html) domain size in the x-direction
    + [`Ly`](/classes/wavevortextransform/ly.html) domain size in the y-direction
    + [`Lz`](/classes/wavevortextransform/lz.html) domain size in the z-direction
    + [`Nj`](/classes/wavevortextransform/nj.html) points in the j-coordinate, `length(z)`
    + [`Nk`](/classes/wavevortextransform/nk.html) points in the k-coordinate, `length(k)`
    + [`Nl`](/classes/wavevortextransform/nl.html) points in the l-coordinate, `length(l)`
    + [`Nx`](/classes/wavevortextransform/nx.html) points in the x-coordinate, `length(x)`
    + [`Ny`](/classes/wavevortextransform/ny.html) points in the y-coordinate, `length(y)`
    + [`Nz`](/classes/wavevortextransform/nz.html) points in the z-coordinate, `length(z)`
    + [`X`](/classes/wavevortextransform/x.html) x-coordinate matrix
    + [`Y`](/classes/wavevortextransform/y.html) y-coordinate matrix
    + [`Z`](/classes/wavevortextransform/z.html) z-coordinate matrix
    + [`j`](/classes/wavevortextransform/j.html) vertical mode number
    + [`k`](/classes/wavevortextransform/k.html) wavenumber-coordinate dimension in the x-direction
    + [`kljGrid`](/classes/wavevortextransform/kljgrid.html) returns the K, L, J coordinate matrices
    + [`l`](/classes/wavevortextransform/l.html) wavenumber-coordinate dimension in the y-direction
    + [`x`](/classes/wavevortextransform/x.html) x-coordinate dimension
    + [`xyzGrid`](/classes/wavevortextransform/xyzgrid.html) returns the X, Y, Z coordinate matrices
    + [`y`](/classes/wavevortextransform/y.html) y-coordinate dimension
    + [`z`](/classes/wavevortextransform/z.html) z-coordinate dimension
  + Stratification
    + [`N0`](/classes/wavevortextransform/n0.html) interior buoyancy frequency at the surface (z=0)
    + [`N2`](/classes/wavevortextransform/n2.html) buoyancy frequency of the mean density
    + [`Nmax`](/classes/wavevortextransform/nmax.html) maximum buoyancy frequency
    + [`dLnN2`](/classes/wavevortextransform/dlnn2.html) d/dz ln N2
    + [`h`](/classes/wavevortextransform/h.html) equivalent depth of each mode
    + [`rho0`](/classes/wavevortextransform/rho0.html) mean density at the surface (z=0)
    + [`rhobar`](/classes/wavevortextransform/rhobar.html) mean density
+ Wave-vortex coefficients
  + [`A0`](/classes/wavevortextransform/a0.html) geostrophic coefficients at reference time t0
  + [`Am`](/classes/wavevortextransform/am.html) negative wave coefficients at reference time t0
  + [`Ap`](/classes/wavevortextransform/ap.html) positive wave coefficients at reference time t0
  + at time $$t$$
    + [`A0t`](/classes/wavevortextransform/a0t.html) geostrophic coefficients at time t
    + [`Amt`](/classes/wavevortextransform/amt.html) negative wave coefficients at time t
    + [`Apt`](/classes/wavevortextransform/apt.html) positive wave coefficients at time t
+ Initial Conditions
  + [`initWithRandomFlow`](/classes/wavevortextransform/initwithrandomflow.html) initialize with a randomized flow
  + [`initWithUVEta`](/classes/wavevortextransform/initwithuveta.html) initialize with fluid variables $$(u,v,\eta)$$
  + [`initWithUVRho`](/classes/wavevortextransform/initwithuvrho.html) initialize with fluid variables $$(u,v,\rho)$$
  + Waves
    + [`addWaveModes`](/classes/wavevortextransform/addwavemodes.html) add amplitudes of the given wave modes
    + [`initWithGMSpectrum`](/classes/wavevortextransform/initwithgmspectrum.html) initialize with a Garrett-Munk spectrum
    + [`initWithSpectralFunction`](/classes/wavevortextransform/initwithspectralfunction.html) initialize the wave spectrum with a given function
    + [`initWithWaveModes`](/classes/wavevortextransform/initwithwavemodes.html) initialize with the given wave modes
    + [`removeAllWaves`](/classes/wavevortextransform/removeallwaves.html) removes all wave from the model, including inertial oscillations
    + [`setWaveModes`](/classes/wavevortextransform/setwavemodes.html) set amplitudes of the given wave modes
    + [`waveCoefficientsFromWaveModes`](/classes/wavevortextransform/wavecoefficientsfromwavemodes.html) Returns the indices (and re-normalized values) of the wave mode appropriate for the Ap, Am matrices.
    + [`waveModesFromWaveCoefficients`](/classes/wavevortextransform/wavemodesfromwavecoefficients.html) Returns normalized amplitudes and phases of all waves
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/wavevortextransform/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/wavevortextransform/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/wavevortextransform/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/wavevortextransform/setinertialmotions.html) set inertial motions
  + Geostrophic Motions
    + [`addGeostrophicStreamfunction`](/classes/wavevortextransform/addgeostrophicstreamfunction.html) add a geostrophic streamfunction to existing geostrophic motions
    + [`initWithGeostrophicStreamfunction`](/classes/wavevortextransform/initwithgeostrophicstreamfunction.html) initialize with a geostrophic streamfunction
    + [`removeAllGeostrophicMotions`](/classes/wavevortextransform/removeallgeostrophicmotions.html) remove all geostrophic motions
    + [`setGeostrophicStreamfunction`](/classes/wavevortextransform/setgeostrophicstreamfunction.html) set a geostrophic streamfunction
+ Energetics
  + [`summarizeEnergyContent`](/classes/wavevortextransform/summarizeenergycontent.html) displays a summary of the energy content of the fluid
  + [`totalEnergy`](/classes/wavevortextransform/totalenergy.html) horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients
  + [`totalEnergySpatiallyIntegrated`](/classes/wavevortextransform/totalenergyspatiallyintegrated.html) horizontally-averaged depth-integrated energy computed in the spatial domain
  + [`totalHydrostaticEnergy`](/classes/wavevortextransform/totalhydrostaticenergy.html) horizontally-averaged depth-integrated energy *without w* computed in the spatial domain
  + Major Constituents
    + [`geostrophicEnergy`](/classes/wavevortextransform/geostrophicenergy.html) total energy, inertial oscillations
    + [`inertialEnergy`](/classes/wavevortextransform/inertialenergy.html) total energy, inertial oscillations
    + [`waveEnergy`](/classes/wavevortextransform/waveenergy.html) total energy, waves
  + Geostrophic Constituents
    + [`geostrophicEnergyBaroclinic`](/classes/wavevortextransform/geostrophicenergybaroclinic.html) total energy, geostrophic, baroclinic
    + [`geostrophicEnergyBarotropic`](/classes/wavevortextransform/geostrophicenergybarotropic.html) total energy, geostrophic, barotropic
  + Inertia-Gravity Wave Constituents
    + [`inertialEnergyBaroclinic`](/classes/wavevortextransform/inertialenergybaroclinic.html) total energy, inertial oscillations, baroclinic
    + [`inertialEnergyBarotropic`](/classes/wavevortextransform/inertialenergybarotropic.html) total energy, inertial oscillations, barotropic
    + [`internalWaveEnergyMinus`](/classes/wavevortextransform/internalwaveenergyminus.html) total energy, internal waves, minus
    + [`internalWaveEnergyPlus`](/classes/wavevortextransform/internalwaveenergyplus.html) total energy, internal waves, positive
  + Multiplicative factors
    + [`A0_HKE_factor`](/classes/wavevortextransform/a0_hke_factor.html) multiplicative factor that multiplies $$A_0$$ to compute horizontal kinetic energy.
    + [`A0_PE_factor`](/classes/wavevortextransform/a0_pe_factor.html) multiplicative factor that multiplies $$A_0$$ to compute potential energy.
    + [`A0_TE_factor`](/classes/wavevortextransform/a0_te_factor.html) multiplicative factor that multiplies $$A_0$$ to compute total energy.
    + [`Apm_TE_factor`](/classes/wavevortextransform/apm_te_factor.html) multiplicative factor that multiplies $$A_\pm$$ to compute total energy.
+ Wave-vortex sorting matrix
  + inverse components ($$S^{-1}$$)
    + [`A0N`](/classes/wavevortextransform/a0n.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_0$$.
    + [`A0U`](/classes/wavevortextransform/a0u.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_0$$.
    + [`A0V`](/classes/wavevortextransform/a0v.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_0$$.
    + [`AmN`](/classes/wavevortextransform/amn.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.
    + [`AmU`](/classes/wavevortextransform/amu.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_m$$.
    + [`AmV`](/classes/wavevortextransform/amv.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_m$$.
    + [`ApN`](/classes/wavevortextransform/apn.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.
    + [`ApU`](/classes/wavevortextransform/apu.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_p$$.
    + [`ApV`](/classes/wavevortextransform/apv.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_p$$.
  + components of $$S$$
    + [`NA0`](/classes/wavevortextransform/na0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{\eta}$$.
    + [`NAm`](/classes/wavevortextransform/nam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{\eta}$$.
    + [`NAp`](/classes/wavevortextransform/nap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{\eta}$$.
    + [`UA0`](/classes/wavevortextransform/ua0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{u}$$.
    + [`UAm`](/classes/wavevortextransform/uam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{u}$$.
    + [`UAp`](/classes/wavevortextransform/uap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{u}$$.
    + [`VA0`](/classes/wavevortextransform/va0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{v}$$.
    + [`VAm`](/classes/wavevortextransform/vam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{v}$$.
    + [`VAp`](/classes/wavevortextransform/vap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{v}$$.
    + [`WAm`](/classes/wavevortextransform/wam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{w}$$.
    + [`WAp`](/classes/wavevortextransform/wap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{w}$$.
+ Other
  + [`EnergeticsByWavenumberAndMode`](/classes/wavevortextransform/energeticsbywavenumberandmode.html) 
  + [`EnergyFluxAtTimeInitial`](/classes/wavevortextransform/energyfluxattimeinitial.html) 
  + [`EnergyFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/energyfluxforflowconstituentsattime.html) 
  + [`EnergyFluxForFlowConstituentsAtTimeInitial`](/classes/wavevortextransform/energyfluxforflowconstituentsattimeinitial.html) 
  + [`ExponentialFilter`](/classes/wavevortextransform/exponentialfilter.html) 
  + [`NonlinearFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/nonlinearfluxforflowconstituentsattime.html) Apply operator T_\omega---defined in (C2) in the manuscript
  + [`WaveVortexTransform`](/classes/wavevortextransform/wavevortextransform.html) These first properties are directly set on initialization
  + [`energyFlux`](/classes/wavevortextransform/energyflux.html) 
  + [`iOmega`](/classes/wavevortextransform/iomega.html) 
  + [`interpolatedFieldAtPositionBadBoundaries`](/classes/wavevortextransform/interpolatedfieldatpositionbadboundaries.html) 
  + [`nonlinearFlux`](/classes/wavevortextransform/nonlinearflux.html) 
  + [`offgridModes`](/classes/wavevortextransform/offgridmodes.html) offgridModes -  subclass should initialize
  + [`ongridModes`](/classes/wavevortextransform/ongridmodes.html) ongridModes -  This is a cached copy
  + [`radialWavenumberAxis`](/classes/wavevortextransform/radialwavenumberaxis.html) Create a reasonable wavenumber axis
  + [`spectralVanishingViscosityFilter`](/classes/wavevortextransform/spectralvanishingviscosityfilter.html) Builds the spectral vanishing viscosity operator
  + [`uMaxGNormRatioForWave`](/classes/wavevortextransform/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`variables`](/classes/wavevortextransform/variables.html) Primary method for accessing the dynamical variables on the
  + [`variablesAtPosition`](/classes/wavevortextransform/variablesatposition.html) Primary method for accessing the dynamical variables on the
  + [`velocityField`](/classes/wavevortextransform/velocityfield.html) Return the velocity field, which is the sum of the gridded
  + [`version`](/classes/wavevortextransform/version.html) 
  + [`waveVortexCoefficientsAtTimeT`](/classes/wavevortextransform/wavevortexcoefficientsattimet.html) 
+ State Variables
  + [`F0`](/classes/wavevortextransform/f0.html) non-linear flux into A0
  + [`Fm`](/classes/wavevortextransform/fm.html) non-linear flux into Am
  + [`Fp`](/classes/wavevortextransform/fp.html) non-linear flux into Ap
  + [`eta`](/classes/wavevortextransform/eta.html) isopycnal deviation
  + [`p`](/classes/wavevortextransform/p.html) pressure anomaly
  + [`qgpv`](/classes/wavevortextransform/qgpv.html) quasigeostrophic potential vorticity
  + [`u`](/classes/wavevortextransform/u.html) x-component of the fluid velocity
  + [`uMax`](/classes/wavevortextransform/umax.html) max horizontal fluid speed
  + [`v`](/classes/wavevortextransform/v.html) y-component of the fluid velocity
  + [`w`](/classes/wavevortextransform/w.html) z-component of the fluid velocity
+ Utility function
  + [`checkHermitian`](/classes/wavevortextransform/checkhermitian.html) Check if the matrix is Hermitian. Report errors.
  + [`extractNonzeroWaveProperties`](/classes/wavevortextransform/extractnonzerowaveproperties.html) Takes a Hermitian matrix and returns the amplitude and phase of nonzero components
  + [`generateHermitianRandomMatrix`](/classes/wavevortextransform/generatehermitianrandommatrix.html) Generate a 3D matrix to be Hermitian, except at k=l=0
  + [`generateRandomFlowState`](/classes/wavevortextransform/generaterandomflowstate.html) Random flow state, separated out by solution type.
  + [`makeHermitian`](/classes/wavevortextransform/makehermitian.html) Forces a 3D matrix to be Hermitian
  + [`nyquistWavenumbers`](/classes/wavevortextransform/nyquistwavenumbers.html) Returns a matrix with 1s at the Nyquist frequencies.
  + [`redundantHermitianCoefficients`](/classes/wavevortextransform/redundanthermitiancoefficients.html) Returns a matrix with 1s at the 'redundant' hermiation indices.
  + Metadata
    + [`addDimensionAnnotations`](/classes/wavevortextransform/adddimensionannotations.html) add one or more WVDimensions
    + [`addOperation`](/classes/wavevortextransform/addoperation.html) add a WVOperation
    + [`addPropertyAnnotations`](/classes/wavevortextransform/addpropertyannotations.html) add a addProperty
    + [`dimensionAnnotationWithName`](/classes/wavevortextransform/dimensionannotationwithname.html) retrieve a WVDimension by name
    + [`operationWithName`](/classes/wavevortextransform/operationwithname.html) retrieve a WVOperation by name
    + [`propertyAnnotationWithName`](/classes/wavevortextransform/propertyannotationwithname.html) retrieve a WVPropertyAnnotation by name
    + [`variableAnnotationWithName`](/classes/wavevortextransform/variableannotationwithname.html) retrieve a WVVariableAnnotation by name
    + [`variableNames`](/classes/wavevortextransform/variablenames.html) retrieve the names of all available variables
+ External (non-gridded) modes
  + [`addExternalWavesWithFrequencies`](/classes/wavevortextransform/addexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given wavenumber
  + [`addExternalWavesWithWavenumbers`](/classes/wavevortextransform/addexternalwaveswithwavenumbers.html) add external (non-gridded) waves with a given wavenumber
  + [`externalVariableFieldsAtTime`](/classes/wavevortextransform/externalvariablefieldsattime.html) Returns the external wave modes at the grid points.
  + [`externalVariablesAtTimePosition`](/classes/wavevortextransform/externalvariablesattimeposition.html) Returns the external wave modes at the grid points.
  + [`fillOutWaveSpectrum`](/classes/wavevortextransform/filloutwavespectrum.html) Add external waves to the model to fill out the spectrum
  + [`removeAllExternalWaves`](/classes/wavevortextransform/removeallexternalwaves.html) remove all external (non-gridded) waves
  + [`setExternalWavesWithFrequencies`](/classes/wavevortextransform/setexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given frequency
  + [`setExternalWavesWithWavenumbers`](/classes/wavevortextransform/setexternalwaveswithwavenumbers.html) set external (non-gridded) waves with a given wavenumber
+ Internal
  + [`addToVariableCache`](/classes/wavevortextransform/addtovariablecache.html) add variable to internal cache, in case it is needed again
  + [`buildTransformationMatrices`](/classes/wavevortextransform/buildtransformationmatrices.html) Part of the internal initialization process where the coefficients for the transformation matrices are constructed.
  + [`clearVariableCache`](/classes/wavevortextransform/clearvariablecache.html) clear the internal cache
  + [`clearVariableCacheOfTimeDependentVariables`](/classes/wavevortextransform/clearvariablecacheoftimedependentvariables.html) clear the internal cache of variables that claim to be time dependent
  + [`defaultDimensionAnnotations`](/classes/wavevortextransform/defaultdimensionannotations.html) return array of TransformDimensions initialized by default
  + [`defaultMethodAnnotations`](/classes/wavevortextransform/defaultmethodannotations.html) return array of WVAnnotations to annotate the methods
  + [`defaultOperations`](/classes/wavevortextransform/defaultoperations.html) return array of WVOperation instances initialized by default
  + [`defaultPropertyAnnotations`](/classes/wavevortextransform/defaultpropertyannotations.html) return array of WVPropertyAnnotation initialized by default
  + [`fetchFromVariableCache`](/classes/wavevortextransform/fetchfromvariablecache.html) retrieve a set of variables from the internal cache
  + [`performTransformOperation`](/classes/wavevortextransform/performtransformoperation.html) computes (runs) the operation
  + [`stateVariables`](/classes/wavevortextransform/statevariables.html) retrieve variables either from cache or by computation
+ Operations
  + Differentiation
    + [`diffX`](/classes/wavevortextransform/diffx.html) differentiate a spatial variable in the x-direction
    + [`diffY`](/classes/wavevortextransform/diffy.html) differentiate a spatial variable in the y-direction
    + [`diffZF`](/classes/wavevortextransform/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes/wavevortextransform/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
  + Transformations
    + [`transformFromSpatialDomainWithF`](/classes/wavevortextransform/transformfromspatialdomainwithf.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the F-modes
    + [`transformFromSpatialDomainWithG`](/classes/wavevortextransform/transformfromspatialdomainwithg.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the G-modes
    + [`transformToRadialWavenumber`](/classes/wavevortextransform/transformtoradialwavenumber.html) transforms in the spectral domain from (k,l,j) to (kRadial,j)
    + [`transformToSpatialDomainWithF`](/classes/wavevortextransform/transformtospatialdomainwithf.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes
    + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithfallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.
    + [`transformToSpatialDomainWithG`](/classes/wavevortextransform/transformtospatialdomainwithg.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes
    + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithgallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.
    + [`transformUVEtaToWaveVortex`](/classes/wavevortextransform/transformuvetatowavevortex.html) transform fluid variables $$(u,v,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.
    + [`transformWaveVortexToUVWEta`](/classes/wavevortextransform/transformwavevortextouvweta.html) transform wave-vortex coefficients $$(A_+,A_-,A_0)$$ to fluid variables $$(u,v,\eta)$$.
+ Masks
  + [`maskForAliasedModes`](/classes/wavevortextransform/maskforaliasedmodes.html) returns a mask with locations of modes that will alias with a quadratic multiplication.
  + [`masksForAllFlowConstituents`](/classes/wavevortextransform/masksforallflowconstituents.html) Returns six 'masks' (matrices with 1s or 0s) indicating where the six
  + [`masksForFlowConstituents`](/classes/wavevortextransform/masksforflowconstituents.html) Returns a sets of 'masks' indicating where different solution types live in the Ap, Am, A0 matrices.
+ Validation and internal unit testing
  + [`validateTransformationMatrices`](/classes/wavevortextransform/validatetransformationmatrices.html) used to confirm if $$S$$ and $$S^{-1}$$ are inverses
+ Write to file
  + [`writeToFile`](/classes/wavevortextransform/writetofile.html) Output the `WaveVortexTransform` to file.


---