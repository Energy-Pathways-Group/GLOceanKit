---
layout: default
title: WVTransform
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 1
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
  + [`waveVortexTransformFromFile`](/classes/wvtransform/wavevortextransformfromfile.html) Initialize a WVTransform instance from an existing file
  + [`waveVortexTransformWithDoubleResolution`](/classes/wvtransform/wavevortextransformwithdoubleresolution.html) create a new WVTransform with double resolution
  + [`waveVortexTransformWithResolution`](/classes/wvtransform/wavevortextransformwithresolution.html) create a new WVTransform with increased resolution
+ Domain attributes
  + [`Omega`](/classes/wvtransform/omega.html) frequency of oscillation of the linear waves
  + [`f`](/classes/wvtransform/f.html) Coriolis parameter
  + [`g`](/classes/wvtransform/g.html) gravity of Earth
  + [`inertialPeriod`](/classes/wvtransform/inertialperiod.html) inertial period
  + [`isBarotropic`](/classes/wvtransform/isbarotropic.html) Boolean indicating whether there is a single (equivalent barotropic) mode
  + [`kRadial`](/classes/wvtransform/kradial.html) isotropic wavenumber dimension
  + [`latitude`](/classes/wvtransform/latitude.html) latitude of the simulation
  + [`t`](/classes/wvtransform/t.html) time of observations
  + [`t0`](/classes/wvtransform/t0.html) reference time of Ap, Am, A0
  + Grid
    + [`J`](/classes/wvtransform/j.html) j-coordinate matrix
    + [`K`](/classes/wvtransform/k.html) k-coordinate matrix
    + [`Kh`](/classes/wvtransform/kh.html) horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$
    + [`L`](/classes/wvtransform/l.html) l-coordinate matrix
    + [`Lx`](/classes/wvtransform/lx.html) domain size in the x-direction
    + [`Ly`](/classes/wvtransform/ly.html) domain size in the y-direction
    + [`Lz`](/classes/wvtransform/lz.html) domain size in the z-direction
    + [`Nj`](/classes/wvtransform/nj.html) points in the j-coordinate, `length(z)`
    + [`Nk`](/classes/wvtransform/nk.html) points in the k-coordinate, `length(k)`
    + [`Nl`](/classes/wvtransform/nl.html) points in the l-coordinate, `length(l)`
    + [`Nx`](/classes/wvtransform/nx.html) points in the x-coordinate, `length(x)`
    + [`Ny`](/classes/wvtransform/ny.html) points in the y-coordinate, `length(y)`
    + [`Nz`](/classes/wvtransform/nz.html) points in the z-coordinate, `length(z)`
    + [`X`](/classes/wvtransform/x.html) x-coordinate matrix
    + [`Y`](/classes/wvtransform/y.html) y-coordinate matrix
    + [`Z`](/classes/wvtransform/z.html) z-coordinate matrix
    + [`j`](/classes/wvtransform/j.html) vertical mode number
    + [`k`](/classes/wvtransform/k.html) wavenumber-coordinate dimension in the x-direction
    + [`kljGrid`](/classes/wvtransform/kljgrid.html) returns the K, L, J coordinate matrices
    + [`l`](/classes/wvtransform/l.html) wavenumber-coordinate dimension in the y-direction
    + [`x`](/classes/wvtransform/x.html) x-coordinate dimension
    + [`xyzGrid`](/classes/wvtransform/xyzgrid.html) returns the X, Y, Z coordinate matrices
    + [`y`](/classes/wvtransform/y.html) y-coordinate dimension
    + [`z`](/classes/wvtransform/z.html) z-coordinate dimension
  + Stratification
    + [`N0`](/classes/wvtransform/n0.html) interior buoyancy frequency at the surface (z=0)
    + [`N2`](/classes/wvtransform/n2.html) buoyancy frequency of the mean density
    + [`Nmax`](/classes/wvtransform/nmax.html) maximum buoyancy frequency
    + [`dLnN2`](/classes/wvtransform/dlnn2.html) d/dz ln N2
    + [`h`](/classes/wvtransform/h.html) equivalent depth of each mode
    + [`rho0`](/classes/wvtransform/rho0.html) mean density at the surface (z=0)
    + [`rhobar`](/classes/wvtransform/rhobar.html) mean density
+ Wave-vortex coefficients
  + [`A0`](/classes/wvtransform/a0.html) geostrophic coefficients at reference time t0
  + [`Am`](/classes/wvtransform/am.html) negative wave coefficients at reference time t0
  + [`Ap`](/classes/wvtransform/ap.html) positive wave coefficients at reference time t0
  + at time $$t$$
    + [`A0t`](/classes/wvtransform/a0t.html) geostrophic coefficients at time t
    + [`Amt`](/classes/wvtransform/amt.html) negative wave coefficients at time t
    + [`Apt`](/classes/wvtransform/apt.html) positive wave coefficients at time t
+ Initial Conditions
  + [`initFromNetCDFFile`](/classes/wvtransform/initfromnetcdffile.html) initialize the flow from a NetCDF file
  + [`initWithRandomFlow`](/classes/wvtransform/initwithrandomflow.html) initialize with a randomized flow
  + [`initWithUVEta`](/classes/wvtransform/initwithuveta.html) initialize with fluid variables $$(u,v,\eta)$$
  + [`initWithUVRho`](/classes/wvtransform/initwithuvrho.html) initialize with fluid variables $$(u,v,\rho)$$
  + [`removeEnergyFromAliasedModes`](/classes/wvtransform/removeenergyfromaliasedmodes.html) remove all energy from aliased modes
  + Waves
    + [`addWaveModes`](/classes/wvtransform/addwavemodes.html) add amplitudes of the given wave modes
    + [`initWithGMSpectrum`](/classes/wvtransform/initwithgmspectrum.html) initialize with a Garrett-Munk spectrum
    + [`initWithSpectralFunction`](/classes/wvtransform/initwithspectralfunction.html) initialize the wave spectrum with a given function
    + [`initWithWaveModes`](/classes/wvtransform/initwithwavemodes.html) initialize with the given wave modes
    + [`removeAllWaves`](/classes/wvtransform/removeallwaves.html) removes all wave from the model, including inertial oscillations
    + [`setWaveModes`](/classes/wvtransform/setwavemodes.html) set amplitudes of the given wave modes
    + [`waveCoefficientsFromWaveModes`](/classes/wvtransform/wavecoefficientsfromwavemodes.html) Returns the indices (and re-normalized values) of the wave mode appropriate for the Ap, Am matrices.
    + [`waveModesFromWaveCoefficients`](/classes/wvtransform/wavemodesfromwavecoefficients.html) Returns normalized amplitudes and phases of all waves
  + Inertial Oscillations
    + [`addInertialMotions`](/classes/wvtransform/addinertialmotions.html) add inertial motions to existing inertial motions
    + [`initWithInertialMotions`](/classes/wvtransform/initwithinertialmotions.html) initialize with inertial motions
    + [`removeAllInertialMotions`](/classes/wvtransform/removeallinertialmotions.html) remove all inertial motions
    + [`setInertialMotions`](/classes/wvtransform/setinertialmotions.html) set inertial motions
  + Geostrophic Motions
    + [`addGeostrophicStreamfunction`](/classes/wvtransform/addgeostrophicstreamfunction.html) add a geostrophic streamfunction to existing geostrophic motions
    + [`initWithGeostrophicStreamfunction`](/classes/wvtransform/initwithgeostrophicstreamfunction.html) initialize with a geostrophic streamfunction
    + [`removeAllGeostrophicMotions`](/classes/wvtransform/removeallgeostrophicmotions.html) remove all geostrophic motions
    + [`setGeostrophicStreamfunction`](/classes/wvtransform/setgeostrophicstreamfunction.html) set a geostrophic streamfunction
+ Energetics
  + [`summarizeEnergyContent`](/classes/wvtransform/summarizeenergycontent.html) displays a summary of the energy content of the fluid
  + [`summarizeModeEnergy`](/classes/wvtransform/summarizemodeenergy.html) List the most energetic modes
  + [`totalEnergy`](/classes/wvtransform/totalenergy.html) horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients
  + [`totalEnergySpatiallyIntegrated`](/classes/wvtransform/totalenergyspatiallyintegrated.html) horizontally-averaged depth-integrated energy computed in the spatial domain
  + [`totalHydrostaticEnergy`](/classes/wvtransform/totalhydrostaticenergy.html) horizontally-averaged depth-integrated energy *without w* computed in the spatial domain
  + Major Constituents
    + [`geostrophicEnergy`](/classes/wvtransform/geostrophicenergy.html) total energy, inertial oscillations
    + [`inertialEnergy`](/classes/wvtransform/inertialenergy.html) total energy, inertial oscillations
    + [`waveEnergy`](/classes/wvtransform/waveenergy.html) total energy, waves
  + Geostrophic Constituents
    + [`geostrophicEnergyBaroclinic`](/classes/wvtransform/geostrophicenergybaroclinic.html) total energy, geostrophic, baroclinic
    + [`geostrophicEnergyBarotropic`](/classes/wvtransform/geostrophicenergybarotropic.html) total energy, geostrophic, barotropic
  + Inertia-Gravity Wave Constituents
    + [`inertialEnergyBaroclinic`](/classes/wvtransform/inertialenergybaroclinic.html) total energy, inertial oscillations, baroclinic
    + [`inertialEnergyBarotropic`](/classes/wvtransform/inertialenergybarotropic.html) total energy, inertial oscillations, barotropic
    + [`internalWaveEnergyMinus`](/classes/wvtransform/internalwaveenergyminus.html) total energy, internal waves, minus
    + [`internalWaveEnergyPlus`](/classes/wvtransform/internalwaveenergyplus.html) total energy, internal waves, positive
  + Multiplicative factors
    + [`A0_HKE_factor`](/classes/wvtransform/a0_hke_factor.html) multiplicative factor that multiplies $$A_0$$ to compute horizontal kinetic energy.
    + [`A0_PE_factor`](/classes/wvtransform/a0_pe_factor.html) multiplicative factor that multiplies $$A_0$$ to compute potential energy.
    + [`A0_TE_factor`](/classes/wvtransform/a0_te_factor.html) multiplicative factor that multiplies $$A_0$$ to compute total energy.
    + [`Apm_TE_factor`](/classes/wvtransform/apm_te_factor.html) multiplicative factor that multiplies $$A_\pm$$ to compute total energy.
+ Wave-vortex sorting matrix
  + inverse components ($$S^{-1}$$)
    + [`A0N`](/classes/wvtransform/a0n.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_0$$.
    + [`A0U`](/classes/wvtransform/a0u.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_0$$.
    + [`A0V`](/classes/wvtransform/a0v.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_0$$.
    + [`AmN`](/classes/wvtransform/amn.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.
    + [`AmU`](/classes/wvtransform/amu.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_m$$.
    + [`AmV`](/classes/wvtransform/amv.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_m$$.
    + [`ApN`](/classes/wvtransform/apn.html) matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.
    + [`ApU`](/classes/wvtransform/apu.html) matrix component that multiplies $$\tilde{u}$$ to compute $$A_p$$.
    + [`ApV`](/classes/wvtransform/apv.html) matrix component that multiplies $$\tilde{v}$$ to compute $$A_p$$.
  + components of $$S$$
    + [`NA0`](/classes/wvtransform/na0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{\eta}$$.
    + [`NAm`](/classes/wvtransform/nam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{\eta}$$.
    + [`NAp`](/classes/wvtransform/nap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{\eta}$$.
    + [`UA0`](/classes/wvtransform/ua0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{u}$$.
    + [`UAm`](/classes/wvtransform/uam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{u}$$.
    + [`UAp`](/classes/wvtransform/uap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{u}$$.
    + [`VA0`](/classes/wvtransform/va0.html) matrix component that multiplies $$A_0$$ to compute $$\tilde{v}$$.
    + [`VAm`](/classes/wvtransform/vam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{v}$$.
    + [`VAp`](/classes/wvtransform/vap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{v}$$.
    + [`WAm`](/classes/wvtransform/wam.html) matrix component that multiplies $$A_m$$ to compute $$\tilde{w}$$.
    + [`WAp`](/classes/wvtransform/wap.html) matrix component that multiplies $$A_p$$ to compute $$\tilde{w}$$.
+ Other
  + [`EnergeticsByWavenumberAndMode`](/classes/wvtransform/energeticsbywavenumberandmode.html) 
  + [`ExponentialFilter`](/classes/wvtransform/exponentialfilter.html) 
  + [`enstrophyFluxFromF0`](/classes/wvtransform/enstrophyfluxfromf0.html) 
  + [`iOmega`](/classes/wvtransform/iomega.html) 
  + [`offgridModes`](/classes/wvtransform/offgridmodes.html) subclass should initialize
  + [`ongridModes`](/classes/wvtransform/ongridmodes.html) This is a cached copy
  + [`qgpvFluxFromF0`](/classes/wvtransform/qgpvfluxfromf0.html) 
  + [`radialWavenumberAxis`](/classes/wvtransform/radialwavenumberaxis.html) Create a reasonable wavenumber axis
  + [`spectralVanishingViscosityFilter`](/classes/wvtransform/spectralvanishingviscosityfilter.html) Builds the spectral vanishing viscosity operator
  + [`uMaxGNormRatioForWave`](/classes/wvtransform/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
  + [`variables`](/classes/wvtransform/variables.html) access the dynamical variables
  + [`variablesAtPosition`](/classes/wvtransform/variablesatposition.html) access the dynamical variables at any position in the domain
  + [`velocityField`](/classes/wvtransform/velocityfield.html) Return the velocity field, which is the sum of the gridded
  + [`version`](/classes/wvtransform/version.html) 
  + [`waveVortexCoefficientsAtTimeT`](/classes/wvtransform/wavevortexcoefficientsattimet.html) 
+ State Variables
  + [`F0`](/classes/wvtransform/f0.html) non-linear flux into A0
  + [`Fm`](/classes/wvtransform/fm.html) non-linear flux into Am
  + [`Fp`](/classes/wvtransform/fp.html) non-linear flux into Ap
  + [`eta`](/classes/wvtransform/eta.html) isopycnal deviation
  + [`p`](/classes/wvtransform/p.html) pressure anomaly
  + [`psi`](/classes/wvtransform/psi.html) geostrophic streamfunction
  + [`qgpv`](/classes/wvtransform/qgpv.html) quasigeostrophic potential vorticity
  + [`rho_prime`](/classes/wvtransform/rho_prime.html) density anomaly
  + [`rho_total`](/classes/wvtransform/rho_total.html) total potential density
  + [`seaSurfaceHeight`](/classes/wvtransform/seasurfaceheight.html) sea-surface height
  + [`seaSurfaceU`](/classes/wvtransform/seasurfaceu.html) x-component of the fluid velocity at the surface
  + [`seaSurfaceV`](/classes/wvtransform/seasurfacev.html) y-component of the fluid velocity at the surface
  + [`u`](/classes/wvtransform/u.html) x-component of the fluid velocity
  + [`uMax`](/classes/wvtransform/umax.html) max horizontal fluid speed
  + [`v`](/classes/wvtransform/v.html) y-component of the fluid velocity
  + [`w`](/classes/wvtransform/w.html) z-component of the fluid velocity
  + [`wMax`](/classes/wvtransform/wmax.html) max vertical fluid speed
+ Internal
  + [`WVTransform`](/classes/wvtransform/wvtransform.html) initialize a WVTransform instance
  + [`addToVariableCache`](/classes/wvtransform/addtovariablecache.html) add variable to internal cache, in case it is needed again
  + [`buildTransformationMatrices`](/classes/wvtransform/buildtransformationmatrices.html) Part of the internal initialization process where the coefficients for the transformation matrices are constructed.
  + [`clearVariableCache`](/classes/wvtransform/clearvariablecache.html) clear the internal cache
  + [`clearVariableCacheOfTimeDependentVariables`](/classes/wvtransform/clearvariablecacheoftimedependentvariables.html) clear the internal cache of variables that claim to be time dependent
  + [`defaultDimensionAnnotations`](/classes/wvtransform/defaultdimensionannotations.html) return array of TransformDimensions initialized by default
  + [`defaultMethodAnnotations`](/classes/wvtransform/defaultmethodannotations.html) return array of WVAnnotations to annotate the methods
  + [`defaultOperations`](/classes/wvtransform/defaultoperations.html) return array of WVOperation instances initialized by default
  + [`defaultPropertyAnnotations`](/classes/wvtransform/defaultpropertyannotations.html) return array of WVPropertyAnnotation initialized by default
  + [`fetchFromVariableCache`](/classes/wvtransform/fetchfromvariablecache.html) retrieve a set of variables from the internal cache
  + [`performOperation`](/classes/wvtransform/performoperation.html) computes (runs) the operation
  + [`performOperationWithName`](/classes/wvtransform/performoperationwithname.html) computes (runs) the operation
  + [`stateVariables`](/classes/wvtransform/statevariables.html) retrieve variables either from cache or by computation
+ Utility function
  + [`checkHermitian`](/classes/wvtransform/checkhermitian.html) Check if the matrix is Hermitian. Report errors.
  + [`extractNonzeroWaveProperties`](/classes/wvtransform/extractnonzerowaveproperties.html) Takes a Hermitian matrix and returns the amplitude and phase of nonzero components
  + [`generateHermitianRandomMatrix`](/classes/wvtransform/generatehermitianrandommatrix.html) Generate a 3D matrix to be Hermitian, except at k=l=0
  + [`generateRandomFlowState`](/classes/wvtransform/generaterandomflowstate.html) Random flow state, separated out by solution type.
  + [`makeHermitian`](/classes/wvtransform/makehermitian.html) Forces a 3D matrix to be Hermitian
  + [`redundantHermitianCoefficients`](/classes/wvtransform/redundanthermitiancoefficients.html) Returns a matrix with 1s at the 'redundant' hermiation indices.
  + Metadata
    + [`addDimensionAnnotations`](/classes/wvtransform/adddimensionannotations.html) add one or more WVDimensions
    + [`addOperation`](/classes/wvtransform/addoperation.html) add a WVOperation
    + [`addPropertyAnnotations`](/classes/wvtransform/addpropertyannotations.html) add a addProperty
    + [`dimensionAnnotationWithName`](/classes/wvtransform/dimensionannotationwithname.html) retrieve a WVDimension by name
    + [`operationWithName`](/classes/wvtransform/operationwithname.html) retrieve a WVOperation by name
    + [`propertyAnnotationWithName`](/classes/wvtransform/propertyannotationwithname.html) retrieve a WVPropertyAnnotation by name
    + [`removeOperation`](/classes/wvtransform/removeoperation.html) remove an existing WVOperation
    + [`variableAnnotationWithName`](/classes/wvtransform/variableannotationwithname.html) retrieve a WVVariableAnnotation by name
    + [`variableNames`](/classes/wvtransform/variablenames.html) retrieve the names of all available variables
+ External (non-gridded) modes
  + [`addExternalWavesWithFrequencies`](/classes/wvtransform/addexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given wavenumber
  + [`addExternalWavesWithWavenumbers`](/classes/wvtransform/addexternalwaveswithwavenumbers.html) add external (non-gridded) waves with a given wavenumber
  + [`externalVariableFieldsAtTime`](/classes/wvtransform/externalvariablefieldsattime.html) Returns the external wave modes at the grid points.
  + [`externalVariablesAtTimePosition`](/classes/wvtransform/externalvariablesattimeposition.html) Returns the external wave modes at the grid points.
  + [`fillOutWaveSpectrum`](/classes/wvtransform/filloutwavespectrum.html) Add external waves to the model to fill out the spectrum
  + [`removeAllExternalWaves`](/classes/wvtransform/removeallexternalwaves.html) remove all external (non-gridded) waves
  + [`setExternalWavesWithFrequencies`](/classes/wvtransform/setexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given frequency
  + [`setExternalWavesWithWavenumbers`](/classes/wvtransform/setexternalwaveswithwavenumbers.html) set external (non-gridded) waves with a given wavenumber
+ Operations
  + Differentiation
    + [`diffX`](/classes/wvtransform/diffx.html) differentiate a spatial variable in the x-direction
    + [`diffY`](/classes/wvtransform/diffy.html) differentiate a spatial variable in the y-direction
    + [`diffZF`](/classes/wvtransform/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
    + [`diffZG`](/classes/wvtransform/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
  + Transformations
    + [`transformFromSpatialDomainWithF`](/classes/wvtransform/transformfromspatialdomainwithf.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the F-modes
    + [`transformFromSpatialDomainWithG`](/classes/wvtransform/transformfromspatialdomainwithg.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the G-modes
    + [`transformToRadialWavenumber`](/classes/wvtransform/transformtoradialwavenumber.html) transforms in the spectral domain from (k,l,j) to (kRadial,j)
    + [`transformToSpatialDomainWithF`](/classes/wvtransform/transformtospatialdomainwithf.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes
    + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wvtransform/transformtospatialdomainwithfallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.
    + [`transformToSpatialDomainWithG`](/classes/wvtransform/transformtospatialdomainwithg.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes
    + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wvtransform/transformtospatialdomainwithgallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.
    + [`transformUVEtaToWaveVortex`](/classes/wvtransform/transformuvetatowavevortex.html) transform fluid variables $$(u,v,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.
    + [`transformWaveVortexToUVWEta`](/classes/wvtransform/transformwavevortextouvweta.html) transform wave-vortex coefficients $$(A_+,A_-,A_0)$$ to fluid variables $$(u,v,\eta)$$.
+ Nonlinear flux and energy transfers
  + [`energyFluxFromNonlinearFlux`](/classes/wvtransform/energyfluxfromnonlinearflux.html) converts nonlinear flux into energy flux
  + [`nonlinearFlux`](/classes/wvtransform/nonlinearflux.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxForFlowConstituents`](/classes/wvtransform/nonlinearfluxforflowconstituents.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxOperation`](/classes/wvtransform/nonlinearfluxoperation.html) The operation responsible for computing the nonlinear flux
  + [`nonlinearFluxWithGradientMasks`](/classes/wvtransform/nonlinearfluxwithgradientmasks.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxWithMask`](/classes/wvtransform/nonlinearfluxwithmask.html) returns the flux of each coefficient as determined by the nonlinear flux
+ Masks
  + [`maskForAliasedModes`](/classes/wvtransform/maskforaliasedmodes.html) returns a mask with locations of modes that will alias with a quadratic multiplication.
  + [`maskForNyquistModes`](/classes/wvtransform/maskfornyquistmodes.html) returns a mask with locations of modes that are not fully resolved
  + [`masksForAllFlowConstituents`](/classes/wvtransform/masksforallflowconstituents.html) Returns six 'masks' (matrices with 1s or 0s) indicating where the six
  + [`masksForFlowConstituents`](/classes/wvtransform/masksforflowconstituents.html) Returns a sets of 'masks' indicating where different solution types live in the Ap, Am, A0 matrices.
+ Validation and internal unit testing
  + [`validateTransformationMatrices`](/classes/wvtransform/validatetransformationmatrices.html) used to confirm if $$S$$ and $$S^{-1}$$ are inverses
+ Write to file
  + [`writeToFile`](/classes/wvtransform/writetofile.html) Output the `WVTransform` to file.


---