---
layout: default
title: WVTransform
parent: WV transform & model
has_children: false
has_toc: false
mathjax: true
---

#  WVTransform

Represents the state of the ocean in terms of energetically orthogonal wave and geostrophic (vortex) solutions


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVTransform < handle</code></pre></div></div>

## Overview
 
 
  The WVTransform subclasses encapsulate data representing the
  state of the ocean at a given instant in time. What makes the
  WVTransform subclasses special is that the state of the ocean
  is represented as energetically independent waves and geostrophic
  motions (vortices). These classes can be queried for any ocean state
  variable including $$u$$, $$v$$, $$w$$, $$\rho$$, $$p$$, but also
  Ertel PV, relative vorticity, or custom defined state variables.
 
  The WVTransform is an abstract class and as such you must
  instatiate one of the concrete subclasses,
 
  + `WVTransformConstantStratification`
  + `WVTransformHydrostatic`
  + `WVTransformSingleMode`
 
                         
  


## Topics
+ Initialization
  + [`waveVortexTransformFromFile`](/classes-transform-and-model/wvtransform/wavevortextransformfromfile.html) Initialize a WVTransform instance from an existing file
  + [`waveVortexTransformWithDoubleResolution`](/classes-transform-and-model/wvtransform/wavevortextransformwithdoubleresolution.html) create a new WVTransform with double resolution
  + [`waveVortexTransformWithResolution`](/classes-transform-and-model/wvtransform/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Domain attributes
  + [`Omega`](/classes-transform-and-model/wvtransform/omega.html) frequency of oscillation of the linear waves
  + [`f`](/classes-transform-and-model/wvtransform/f.html) Coriolis parameter
  + [`g`](/classes-transform-and-model/wvtransform/g.html) gravity of Earth
  + [`inertialPeriod`](/classes-transform-and-model/wvtransform/inertialperiod.html) inertial period
  + [`isBarotropic`](/classes-transform-and-model/wvtransform/isbarotropic.html) Boolean indicating whether there is a single (equivalent barotropic) mode
  + [`kRadial`](/classes-transform-and-model/wvtransform/kradial.html) isotropic wavenumber dimension
  + [`latitude`](/classes-transform-and-model/wvtransform/latitude.html) latitude of the simulation
  + [`t`](/classes-transform-and-model/wvtransform/t.html) time of observations
  + [`t0`](/classes-transform-and-model/wvtransform/t0.html) reference time of Ap, Am, A0
  + Grid
    + [`J`](/classes-transform-and-model/wvtransform/j.html) j-coordinate matrix
    + [`K`](/classes-transform-and-model/wvtransform/k.html) k-coordinate matrix
    + [`Kh`](/classes-transform-and-model/wvtransform/kh.html) horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$
    + [`L`](/classes-transform-and-model/wvtransform/l.html) l-coordinate matrix
    + [`Lx`](/classes-transform-and-model/wvtransform/lx.html) domain size in the x-direction
    + [`Ly`](/classes-transform-and-model/wvtransform/ly.html) domain size in the y-direction
    + [`Lz`](/classes-transform-and-model/wvtransform/lz.html) domain size in the z-direction
    + [`Nj`](/classes-transform-and-model/wvtransform/nj.html) points in the j-coordinate, `length(z)`
    + [`Nk`](/classes-transform-and-model/wvtransform/nk.html) points in the k-coordinate, `length(k)`
    + [`Nl`](/classes-transform-and-model/wvtransform/nl.html) points in the l-coordinate, `length(l)`
    + [`Nx`](/classes-transform-and-model/wvtransform/nx.html) points in the x-coordinate, `length(x)`
    + [`Ny`](/classes-transform-and-model/wvtransform/ny.html) points in the y-coordinate, `length(y)`
    + [`Nz`](/classes-transform-and-model/wvtransform/nz.html) points in the z-coordinate, `length(z)`
    + [`X`](/classes-transform-and-model/wvtransform/x.html) x-coordinate matrix
    + [`Y`](/classes-transform-and-model/wvtransform/y.html) y-coordinate matrix
    + [`Z`](/classes-transform-and-model/wvtransform/z.html) z-coordinate matrix
    + [`j`](/classes-transform-and-model/wvtransform/j.html) vertical mode number
    + [`k`](/classes-transform-and-model/wvtransform/k.html) wavenumber-coordinate dimension in the x-direction
    + [`kljGrid`](/classes-transform-and-model/wvtransform/kljgrid.html) returns the K, L, J coordinate matrices
    + [`l`](/classes-transform-and-model/wvtransform/l.html) wavenumber-coordinate dimension in the y-direction
    + [`x`](/classes-transform-and-model/wvtransform/x.html) x-coordinate dimension
    + [`xyzGrid`](/classes-transform-and-model/wvtransform/xyzgrid.html) returns the X, Y, Z coordinate matrices
    + [`y`](/classes-transform-and-model/wvtransform/y.html) y-coordinate dimension
    + [`z`](/classes-transform-and-model/wvtransform/z.html) z-coordinate dimension
  + Stratification
    + [`N0`](/classes-transform-and-model/wvtransform/n0.html) interior buoyancy frequency at the surface (z=0)
    + [`N2`](/classes-transform-and-model/wvtransform/n2.html) buoyancy frequency of the mean density
    + [`Nmax`](/classes-transform-and-model/wvtransform/nmax.html) maximum buoyancy frequency
    + [`dLnN2`](/classes-transform-and-model/wvtransform/dlnn2.html) d/dz ln N2
    + [`h`](/classes-transform-and-model/wvtransform/h.html) equivalent depth of each mode
    + [`rho0`](/classes-transform-and-model/wvtransform/rho0.html) mean density at the surface (z=0)
    + [`rhobar`](/classes-transform-and-model/wvtransform/rhobar.html) mean density
+ Wave-vortex coefficients
  + [`A0`](/classes-transform-and-model/wvtransform/a0.html) geostrophic coefficients at reference time t0
  + [`Am`](/classes-transform-and-model/wvtransform/am.html) negative wave coefficients at reference time t0
  + [`Ap`](/classes-transform-and-model/wvtransform/ap.html) positive wave coefficients at reference time t0
  + at time $$t$$
    + [`A0t`](/classes-transform-and-model/wvtransform/a0t.html) geostrophic coefficients at time t
    + [`Amt`](/classes-transform-and-model/wvtransform/amt.html) negative wave coefficients at time t
    + [`Apt`](/classes-transform-and-model/wvtransform/apt.html) positive wave coefficients at time t
+ Initial Conditions
  + [`initFromNetCDFFile`](/classes-transform-and-model/wvtransform/initfromnetcdffile.html) initialize the flow from a NetCDF file
  + [`initWithRandomFlow`](/classes-transform-and-model/wvtransform/initwithrandomflow.html) initialize with a randomized flow
  + [`initWithUVEta`](/classes-transform-and-model/wvtransform/initwithuveta.html) initialize with fluid variables $$(u,v,\eta)$$
  + [`initWithUVRho`](/classes-transform-and-model/wvtransform/initwithuvrho.html) initialize with fluid variables $$(u,v,\rho)$$
  + [`removeEnergyFromAliasedModes`](/classes-transform-and-model/wvtransform/removeenergyfromaliasedmodes.html) remove all energy from aliased modes
  + Waves
    + [`addWaveModes`](/classes-transform-and-model/wvtransform/addwavemodes.html) add amplitudes of the given wave modes
    + [`initWithGMSpectrum`](/classes-transform-and-model/wvtransform/initwithgmspectrum.html) initialize with a Garrett-Munk spectrum
    + [`initWithSpectralFunction`](/classes-transform-and-model/wvtransform/initwithspectralfunction.html) initialize the wave spectrum with a given function
    + [`initWithWaveModes`](/classes-transform-and-model/wvtransform/initwithwavemodes.html) initialize with the given wave modes
    + [`removeAllWaves`](/classes-transform-and-model/wvtransform/removeallwaves.html) removes all wave from the model, including inertial oscillations
    + [`setWaveModes`](/classes-transform-and-model/wvtransform/setwavemodes.html) set amplitudes of the given wave modes
    + [`waveCoefficientsFromWaveModes`](/classes-transform-and-model/wvtransform/wavecoefficientsfromwavemodes.html) Returns the indices (and re-normalized values) of the wave mode appropriate for the Ap, Am matrices.
    + [`waveModesFromWaveCoefficients`](/classes-transform-and-model/wvtransform/wavemodesfromwavecoefficients.html) Returns normalized amplitudes and phases of all waves
  + Inertial Oscillations
    + [`addInertialMotions`](/classes-transform-and-model/wvtransform/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes-transform-and-model/wvtransform/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes-transform-and-model/wvtransform/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes-transform-and-model/wvtransform/setinertialmotions.html) set inertial motions
  + Geostrophic Motions
    + [`addGeostrophicStreamfunction`](/classes-transform-and-model/wvtransform/addgeostrophicstreamfunction.html) add a geostrophic streamfunction to existing geostrophic motions
    + [`initWithGeostrophicStreamfunction`](/classes-transform-and-model/wvtransform/initwithgeostrophicstreamfunction.html) initialize with a geostrophic streamfunction
    + [`removeAllGeostrophicMotions`](/classes-transform-and-model/wvtransform/removeallgeostrophicmotions.html) remove all geostrophic motions
    + [`setGeostrophicStreamfunction`](/classes-transform-and-model/wvtransform/setgeostrophicstreamfunction.html) set a geostrophic streamfunction
+ Energetics
  + [`summarizeEnergyContent`](/classes-transform-and-model/wvtransform/summarizeenergycontent.html) displays a summary of the energy content of the fluid
  + [`summarizeModeEnergy`](/classes-transform-and-model/wvtransform/summarizemodeenergy.html) List the most energetic modes
  + [`totalEnergy`](/classes-transform-and-model/wvtransform/totalenergy.html) horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients
  + [`totalEnergySpatiallyIntegrated`](/classes-transform-and-model/wvtransform/totalenergyspatiallyintegrated.html) horizontally-averaged depth-integrated energy computed in the spatial domain
  + [`totalHydrostaticEnergy`](/classes-transform-and-model/wvtransform/totalhydrostaticenergy.html) horizontally-averaged depth-integrated energy *without w* computed in the spatial domain
  + Major Constituents
    + [`geostrophicEnergy`](/classes-transform-and-model/wvtransform/geostrophicenergy.html) total energy, inertial oscillations
    + [`inertialEnergy`](/classes-transform-and-model/wvtransform/inertialenergy.html) total energy, inertial oscillations
    + [`waveEnergy`](/classes-transform-and-model/wvtransform/waveenergy.html) total energy, waves
  + Geostrophic Constituents
    + [`geostrophicEnergyBaroclinic`](/classes-transform-and-model/wvtransform/geostrophicenergybaroclinic.html) total energy, geostrophic, baroclinic
    + [`geostrophicEnergyBarotropic`](/classes-transform-and-model/wvtransform/geostrophicenergybarotropic.html) total energy, geostrophic, barotropic
  + Inertia-Gravity Wave Constituents
    + [`inertialEnergyBaroclinic`](/classes-transform-and-model/wvtransform/inertialenergybaroclinic.html) total energy, inertial oscillations, baroclinic
    + [`inertialEnergyBarotropic`](/classes-transform-and-model/wvtransform/inertialenergybarotropic.html) total energy, inertial oscillations, barotropic
    + [`internalWaveEnergyMinus`](/classes-transform-and-model/wvtransform/internalwaveenergyminus.html) total energy, internal waves, minus
    + [`internalWaveEnergyPlus`](/classes-transform-and-model/wvtransform/internalwaveenergyplus.html) total energy, internal waves, positive
  + Multiplicative factors
    + [`A0_HKE_factor`](/classes-transform-and-model/wvtransform/a0_hke_factor.html) multiplicative factor that multiplies $$A_0$$ to compute horizontal kinetic energy.
    + [`A0_PE_factor`](/classes-transform-and-model/wvtransform/a0_pe_factor.html) multiplicative factor that multiplies $$A_0$$ to compute potential energy.
    + [`A0_TE_factor`](/classes-transform-and-model/wvtransform/a0_te_factor.html) multiplicative factor that multiplies $$A_0$$ to compute total energy.
    + [`Apm_TE_factor`](/classes-transform-and-model/wvtransform/apm_te_factor.html) multiplicative factor that multiplies $$A_\pm$$ to compute total energy.
+ Wave-vortex sorting matrix
  + inverse components ($$S^{-1}$$)
    + [`A0N`](/classes-transform-and-model/wvtransform/a0n.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_0$$.
    + [`A0U`](/classes-transform-and-model/wvtransform/a0u.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_0$$.
    + [`A0V`](/classes-transform-and-model/wvtransform/a0v.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_0$$.
    + [`AmN`](/classes-transform-and-model/wvtransform/amn.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.
    + [`AmU`](/classes-transform-and-model/wvtransform/amu.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_m$$.
    + [`AmV`](/classes-transform-and-model/wvtransform/amv.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_m$$.
    + [`ApN`](/classes-transform-and-model/wvtransform/apn.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.
    + [`ApU`](/classes-transform-and-model/wvtransform/apu.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_p$$.
    + [`ApV`](/classes-transform-and-model/wvtransform/apv.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_p$$.
  + components of $$S$$
    + [`NA0`](/classes-transform-and-model/wvtransform/na0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{\eta}$$.
    + [`NAm`](/classes-transform-and-model/wvtransform/nam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{\eta}$$.
    + [`NAp`](/classes-transform-and-model/wvtransform/nap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{\eta}$$.
    + [`UA0`](/classes-transform-and-model/wvtransform/ua0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{u}$$.
    + [`UAm`](/classes-transform-and-model/wvtransform/uam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{u}$$.
    + [`UAp`](/classes-transform-and-model/wvtransform/uap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{u}$$.
    + [`VA0`](/classes-transform-and-model/wvtransform/va0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{v}$$.
    + [`VAm`](/classes-transform-and-model/wvtransform/vam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{v}$$.
    + [`VAp`](/classes-transform-and-model/wvtransform/vap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{v}$$.
    + [`WAm`](/classes-transform-and-model/wvtransform/wam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{w}$$.
    + [`WAp`](/classes-transform-and-model/wvtransform/wap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{w}$$.
+ Other
  + [`EnergeticsByWavenumberAndMode`](/classes-transform-and-model/wvtransform/energeticsbywavenumberandmode.html) 
  + [`ExponentialFilter`](/classes-transform-and-model/wvtransform/exponentialfilter.html) 
  + [`enstrophyFluxFromF0`](/classes-transform-and-model/wvtransform/enstrophyfluxfromf0.html) 
  + [`iOmega`](/classes-transform-and-model/wvtransform/iomega.html) 
  + [`offgridModes`](/classes-transform-and-model/wvtransform/offgridmodes.html) subclass should initialize
  + [`ongridModes`](/classes-transform-and-model/wvtransform/ongridmodes.html) This is a cached copy
  + [`qgpvFluxFromF0`](/classes-transform-and-model/wvtransform/qgpvfluxfromf0.html) 
  + [`radialWavenumberAxis`](/classes-transform-and-model/wvtransform/radialwavenumberaxis.html) Create a reasonable wavenumber axis
  + [`spectralVanishingViscosityFilter`](/classes-transform-and-model/wvtransform/spectralvanishingviscosityfilter.html) Builds the spectral vanishing viscosity operator
  + [`uMaxGNormRatioForWave`](/classes-transform-and-model/wvtransform/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`variables`](/classes-transform-and-model/wvtransform/variables.html) access the dynamical variables
  + [`variablesAtPosition`](/classes-transform-and-model/wvtransform/variablesatposition.html) access the dynamical variables at any position in the domain
  + [`velocityField`](/classes-transform-and-model/wvtransform/velocityfield.html) Return the velocity field, which is the sum of the gridded
  + [`version`](/classes-transform-and-model/wvtransform/version.html) 
  + [`waveVortexCoefficientsAtTimeT`](/classes-transform-and-model/wvtransform/wavevortexcoefficientsattimet.html) 
+ State Variables
  + [`F0`](/classes-transform-and-model/wvtransform/f0.html) non-linear flux into A0
  + [`Fm`](/classes-transform-and-model/wvtransform/fm.html) non-linear flux into Am
  + [`Fp`](/classes-transform-and-model/wvtransform/fp.html) non-linear flux into Ap
  + [`eta`](/classes-transform-and-model/wvtransform/eta.html) isopycnal deviation
  + [`p`](/classes-transform-and-model/wvtransform/p.html) pressure anomaly
  + [`psi`](/classes-transform-and-model/wvtransform/psi.html) geostrophic streamfunction
  + [`qgpv`](/classes-transform-and-model/wvtransform/qgpv.html) quasigeostrophic potential vorticity
  + [`rho_prime`](/classes-transform-and-model/wvtransform/rho_prime.html) density anomaly
  + [`rho_total`](/classes-transform-and-model/wvtransform/rho_total.html) total potential density
  + [`seaSurfaceHeight`](/classes-transform-and-model/wvtransform/seasurfaceheight.html) sea-surface height
  + [`seaSurfaceU`](/classes-transform-and-model/wvtransform/seasurfaceu.html) x-component of the fluid velocity at the surface
  + [`seaSurfaceV`](/classes-transform-and-model/wvtransform/seasurfacev.html) y-component of the fluid velocity at the surface
  + [`u`](/classes-transform-and-model/wvtransform/u.html) x-component of the fluid velocity
  + [`uMax`](/classes-transform-and-model/wvtransform/umax.html) max horizontal fluid speed
  + [`v`](/classes-transform-and-model/wvtransform/v.html) y-component of the fluid velocity
  + [`w`](/classes-transform-and-model/wvtransform/w.html) z-component of the fluid velocity
  + [`wMax`](/classes-transform-and-model/wvtransform/wmax.html) max vertical fluid speed
+ Internal
  + [`WVTransform`](/classes-transform-and-model/wvtransform/wvtransform.html) initialize a WVTransform instance
  + [`addToVariableCache`](/classes-transform-and-model/wvtransform/addtovariablecache.html) add variable to internal cache, in case it is needed again
  + [`buildTransformationMatrices`](/classes-transform-and-model/wvtransform/buildtransformationmatrices.html) Part of the internal initialization process where the coefficients for the transformation matrices are constructed.
  + [`clearVariableCache`](/classes-transform-and-model/wvtransform/clearvariablecache.html) clear the internal cache
  + [`clearVariableCacheOfTimeDependentVariables`](/classes-transform-and-model/wvtransform/clearvariablecacheoftimedependentvariables.html) clear the internal cache of variables that claim to be time dependent
  + [`defaultDimensionAnnotations`](/classes-transform-and-model/wvtransform/defaultdimensionannotations.html) return array of TransformDimensions initialized by default
  + [`defaultMethodAnnotations`](/classes-transform-and-model/wvtransform/defaultmethodannotations.html) return array of WVAnnotations to annotate the methods
  + [`defaultOperations`](/classes-transform-and-model/wvtransform/defaultoperations.html) return array of WVOperation instances initialized by default
  + [`defaultPropertyAnnotations`](/classes-transform-and-model/wvtransform/defaultpropertyannotations.html) return array of WVPropertyAnnotation initialized by default
  + [`fetchFromVariableCache`](/classes-transform-and-model/wvtransform/fetchfromvariablecache.html) retrieve a set of variables from the internal cache
  + [`performOperation`](/classes-transform-and-model/wvtransform/performoperation.html) computes (runs) the operation
  + [`performOperationWithName`](/classes-transform-and-model/wvtransform/performoperationwithname.html) computes (runs) the operation
  + [`stateVariables`](/classes-transform-and-model/wvtransform/statevariables.html) retrieve variables either from cache or by computation
+ Utility function
  + [`checkHermitian`](/classes-transform-and-model/wvtransform/checkhermitian.html) Check if the matrix is Hermitian. Report errors.
  + [`extractNonzeroWaveProperties`](/classes-transform-and-model/wvtransform/extractnonzerowaveproperties.html) Takes a Hermitian matrix and returns the amplitude and phase of nonzero components
  + [`generateHermitianRandomMatrix`](/classes-transform-and-model/wvtransform/generatehermitianrandommatrix.html) Generate a 3D matrix to be Hermitian, except at k=l=0
  + [`generateRandomFlowState`](/classes-transform-and-model/wvtransform/generaterandomflowstate.html) Random flow state, separated out by solution type.
  + [`makeHermitian`](/classes-transform-and-model/wvtransform/makehermitian.html) Forces a 3D matrix to be Hermitian
  + [`redundantHermitianCoefficients`](/classes-transform-and-model/wvtransform/redundanthermitiancoefficients.html) Returns a matrix with 1s at the 'redundant' hermiation indices.
  + Metadata
    + [`addDimensionAnnotations`](/classes-transform-and-model/wvtransform/adddimensionannotations.html) add one or more WVDimensions
    + [`addOperation`](/classes-transform-and-model/wvtransform/addoperation.html) add a WVOperation
    + [`addPropertyAnnotations`](/classes-transform-and-model/wvtransform/addpropertyannotations.html) add a addProperty
    + [`dimensionAnnotationWithName`](/classes-transform-and-model/wvtransform/dimensionannotationwithname.html) retrieve a WVDimension by name
    + [`operationWithName`](/classes-transform-and-model/wvtransform/operationwithname.html) retrieve a WVOperation by name
    + [`propertyAnnotationWithName`](/classes-transform-and-model/wvtransform/propertyannotationwithname.html) retrieve a WVPropertyAnnotation by name
    + [`removeOperation`](/classes-transform-and-model/wvtransform/removeoperation.html) remove an existing WVOperation
    + [`variableAnnotationWithName`](/classes-transform-and-model/wvtransform/variableannotationwithname.html) retrieve a WVVariableAnnotation by name
    + [`variableNames`](/classes-transform-and-model/wvtransform/variablenames.html) retrieve the names of all available variables
+ External (non-gridded) modes
  + [`addExternalWavesWithFrequencies`](/classes-transform-and-model/wvtransform/addexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given wavenumber
  + [`addExternalWavesWithWavenumbers`](/classes-transform-and-model/wvtransform/addexternalwaveswithwavenumbers.html) add external (non-gridded) waves with a given wavenumber
  + [`externalVariableFieldsAtTime`](/classes-transform-and-model/wvtransform/externalvariablefieldsattime.html) Returns the external wave modes at the grid points.
  + [`externalVariablesAtTimePosition`](/classes-transform-and-model/wvtransform/externalvariablesattimeposition.html) Returns the external wave modes at the grid points.
  + [`fillOutWaveSpectrum`](/classes-transform-and-model/wvtransform/filloutwavespectrum.html) Add external waves to the model to fill out the spectrum
  + [`removeAllExternalWaves`](/classes-transform-and-model/wvtransform/removeallexternalwaves.html) remove all external (non-gridded) waves
  + [`setExternalWavesWithFrequencies`](/classes-transform-and-model/wvtransform/setexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given frequency
  + [`setExternalWavesWithWavenumbers`](/classes-transform-and-model/wvtransform/setexternalwaveswithwavenumbers.html) set external (non-gridded) waves with a given wavenumber
+ Operations
  + Differentiation
    + [`diffX`](/classes-transform-and-model/wvtransform/diffx.html) differentiate a spatial variable in the x-direction
    + [`diffY`](/classes-transform-and-model/wvtransform/diffy.html) differentiate a spatial variable in the y-direction
    + [`diffZF`](/classes-transform-and-model/wvtransform/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes-transform-and-model/wvtransform/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
  + Transformations
    + [`transformFromSpatialDomainWithF`](/classes-transform-and-model/wvtransform/transformfromspatialdomainwithf.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the F-modes
    + [`transformFromSpatialDomainWithG`](/classes-transform-and-model/wvtransform/transformfromspatialdomainwithg.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the G-modes
    + [`transformToRadialWavenumber`](/classes-transform-and-model/wvtransform/transformtoradialwavenumber.html) transforms in the spectral domain from (k,l,j) to (kRadial,j)
    + [`transformToSpatialDomainWithF`](/classes-transform-and-model/wvtransform/transformtospatialdomainwithf.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes
    + [`transformToSpatialDomainWithFAllDerivatives`](/classes-transform-and-model/wvtransform/transformtospatialdomainwithfallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.
    + [`transformToSpatialDomainWithG`](/classes-transform-and-model/wvtransform/transformtospatialdomainwithg.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes
    + [`transformToSpatialDomainWithGAllDerivatives`](/classes-transform-and-model/wvtransform/transformtospatialdomainwithgallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.
    + [`transformUVEtaToWaveVortex`](/classes-transform-and-model/wvtransform/transformuvetatowavevortex.html) transform fluid variables $$(u,v,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.
    + [`transformWaveVortexToUVWEta`](/classes-transform-and-model/wvtransform/transformwavevortextouvweta.html) transform wave-vortex coefficients $$(A_+,A_-,A_0)$$ to fluid variables $$(u,v,\eta)$$.
+ Nonlinear flux and energy transfers
  + [`energyFluxFromNonlinearFlux`](/classes-transform-and-model/wvtransform/energyfluxfromnonlinearflux.html) converts nonlinear flux into energy flux
  + [`nonlinearFlux`](/classes-transform-and-model/wvtransform/nonlinearflux.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxForFlowConstituents`](/classes-transform-and-model/wvtransform/nonlinearfluxforflowconstituents.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxOperation`](/classes-transform-and-model/wvtransform/nonlinearfluxoperation.html) The operation responsible for computing the nonlinear flux
  + [`nonlinearFluxWithGradientMasks`](/classes-transform-and-model/wvtransform/nonlinearfluxwithgradientmasks.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxWithMask`](/classes-transform-and-model/wvtransform/nonlinearfluxwithmask.html) returns the flux of each coefficient as determined by the nonlinear flux
+ Masks
  + [`maskForAliasedModes`](/classes-transform-and-model/wvtransform/maskforaliasedmodes.html) returns a mask with locations of modes that will alias with a quadratic multiplication.
  + [`maskForNyquistModes`](/classes-transform-and-model/wvtransform/maskfornyquistmodes.html) returns a mask with locations of modes that are not fully resolved
  + [`masksForAllFlowConstituents`](/classes-transform-and-model/wvtransform/masksforallflowconstituents.html) Returns six 'masks' (matrices with 1s or 0s) indicating where the six
  + [`masksForFlowConstituents`](/classes-transform-and-model/wvtransform/masksforflowconstituents.html) Returns a sets of 'masks' indicating where different solution types live in the Ap, Am, A0 matrices.
+ Validation and internal unit testing
  + [`validateTransformationMatrices`](/classes-transform-and-model/wvtransform/validatetransformationmatrices.html) used to confirm if $$S$$ and $$S^{-1}$$ are inverses
+ Write to file
  + [`writeToFile`](/classes-transform-and-model/wvtransform/writetofile.html) Output the `WVTransform` to file.


---