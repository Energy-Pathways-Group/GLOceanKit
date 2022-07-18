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

## Discussion
 
 
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
 
                        


## Topics
+ Initialization
  + Other
    + [`transformFromFile`](/classes/wavevortextransform/transformfromfile.html) Initialize a WaveVortexTransform instance from an existing file
+ Domain attributes
  + Grid
    + [`Kh`](/classes/wavevortextransform/kh.html) horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$
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
    + [`l`](/classes/wavevortextransform/l.html) wavenumber-coordinate dimension in the y-direction
    + [`x`](/classes/wavevortextransform/x.html) x-coordinate dimension
    + [`y`](/classes/wavevortextransform/y.html) y-coordinate dimension
    + [`z`](/classes/wavevortextransform/z.html) z-coordinate dimension
  + Other
    + [`F0`](/classes/wavevortextransform/f0.html) non-linear flux into A0
    + [`Omega`](/classes/wavevortextransform/omega.html) frequency of oscillation of the linear waves
    + [`f0`](/classes/wavevortextransform/f0.html) Coriolis parameter
    + [`inertialPeriod`](/classes/wavevortextransform/inertialperiod.html) inertial period
    + [`isBarotropic`](/classes/wavevortextransform/isbarotropic.html) Boolean indicating whether there is a single (equivalent barotropic) mode
    + [`kRadial`](/classes/wavevortextransform/kradial.html) isotropic wavenumber dimension
    + [`latitude`](/classes/wavevortextransform/latitude.html) latitude of the simulation
    + [`t`](/classes/wavevortextransform/t.html) time of observations
    + [`t0`](/classes/wavevortextransform/t0.html) reference time of Ap, Am, A0
  + Stratification
    + [`N0`](/classes/wavevortextransform/n0.html) interior buoyancy frequency at the surface (z=0)
    + [`N2`](/classes/wavevortextransform/n2.html) buoyancy frequency of the mean density
    + [`Nmax`](/classes/wavevortextransform/nmax.html) maximum buoyancy frequency
    + [`dLnN2`](/classes/wavevortextransform/dlnn2.html) d/dz ln N2
    + [`h`](/classes/wavevortextransform/h.html) equivalent depth of each mode
    + [`rho0`](/classes/wavevortextransform/rho0.html) mean density at the surface (z=0)
    + [`rhobar`](/classes/wavevortextransform/rhobar.html) mean density
+ Wave-vortex coefficients
  + Other
    + [`A0`](/classes/wavevortextransform/a0.html) geostrophic coefficients at reference time t0
    + [`Am`](/classes/wavevortextransform/am.html) negative wave coefficients at reference time t0
    + [`Ap`](/classes/wavevortextransform/ap.html) positive wave coefficients at reference time t0
  + at time $$t$$
    + [`A0t`](/classes/wavevortextransform/a0t.html) geostrophic coefficients at time t
    + [`Amt`](/classes/wavevortextransform/amt.html) negative wave coefficients at time t
    + [`Apt`](/classes/wavevortextransform/apt.html) positive wave coefficients at time t
+ Initial Conditions
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
  + Other
    + [`totalEnergy`](/classes/wavevortextransform/totalenergy.html) horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients
    + [`totalEnergySpatiallyIntegrated`](/classes/wavevortextransform/totalenergyspatiallyintegrated.html) horizontally-averaged depth-integrated energy computed in the spatial domain
    + [`totalHydrostaticEnergy`](/classes/wavevortextransform/totalhydrostaticenergy.html) horizontally-averaged depth-integrated energy *without w* computed in the spatial domain
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
  + Other
    + [`EnergeticsByWavenumberAndMode`](/classes/wavevortextransform/energeticsbywavenumberandmode.html) 
    + [`EnergyFluxAtTimeInitial`](/classes/wavevortextransform/energyfluxattimeinitial.html) 
    + [`EnergyFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/energyfluxforflowconstituentsattime.html) 
    + [`EnergyFluxForFlowConstituentsAtTimeInitial`](/classes/wavevortextransform/energyfluxforflowconstituentsattimeinitial.html) 
    + [`ExponentialFilter`](/classes/wavevortextransform/exponentialfilter.html) 
    + [`InitializeWithRandomFlowState`](/classes/wavevortextransform/initializewithrandomflowstate.html) 
    + [`NonlinearFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/nonlinearfluxforflowconstituentsattime.html) Apply operator T_\omega---defined in (C2) in the manuscript
    + [`WaveVortexTransform`](/classes/wavevortextransform/wavevortextransform.html) These first properties are directly set on initialization
    + [`buildTransformationMatrices`](/classes/wavevortextransform/buildtransformationmatrices.html) Build wavenumbers
    + [`diffZF`](/classes/wavevortextransform/diffzf.html) 
    + [`diffZG`](/classes/wavevortextransform/diffzg.html) 
    + [`energyFlux`](/classes/wavevortextransform/energyflux.html) 
    + [`g`](/classes/wavevortextransform/g.html) gravity of Earth
    + [`generateRandomFlowState`](/classes/wavevortextransform/generaterandomflowstate.html) Random flow state, separated out by solution type.
    + [`grid`](/classes/wavevortextransform/grid.html) 
    + [`iOmega`](/classes/wavevortextransform/iomega.html) 
    + [`interpolatedFieldAtPositionBadBoundaries`](/classes/wavevortextransform/interpolatedfieldatpositionbadboundaries.html) 
    + [`nonlinearFlux`](/classes/wavevortextransform/nonlinearflux.html) 
    + [`offgridModes`](/classes/wavevortextransform/offgridmodes.html) offgridModes -  subclass should initialize
    + [`ongridModes`](/classes/wavevortextransform/ongridmodes.html) ongridModes -  This is a cached copy
    + [`radialWavenumberAxis`](/classes/wavevortextransform/radialwavenumberaxis.html) Create a reasonable wavenumber axis
    + [`rebuildTransformationMatrices`](/classes/wavevortextransform/rebuildtransformationmatrices.html) 
    + [`spectralVanishingViscosityFilter`](/classes/wavevortextransform/spectralvanishingviscosityfilter.html) Builds the spectral vanishing viscosity operator
    + [`stateVariableWithName`](/classes/wavevortextransform/statevariablewithname.html) 
    + [`summarizeEnergyContent`](/classes/wavevortextransform/summarizeenergycontent.html) 
    + [`timeDependentStateVariables`](/classes/wavevortextransform/timedependentstatevariables.html) 
    + [`transformDimensionWithName`](/classes/wavevortextransform/transformdimensionwithname.html) 
    + [`transformOperationWithName`](/classes/wavevortextransform/transformoperationwithname.html) 
    + [`transformPropertyWithName`](/classes/wavevortextransform/transformpropertywithname.html) 
    + [`transformToRadialWavenumber`](/classes/wavevortextransform/transformtoradialwavenumber.html) 
    + [`transformUVEtaToWaveVortex`](/classes/wavevortextransform/transformuvetatowavevortex.html) This is the 'S^{-1}' operator (C5) in the manuscript
    + [`transformWaveVortexToUVWEta`](/classes/wavevortextransform/transformwavevortextouvweta.html) 
    + [`transformWithDoubleResolution`](/classes/wavevortextransform/transformwithdoubleresolution.html) 
    + [`transformWithResolution`](/classes/wavevortextransform/transformwithresolution.html) 
    + [`uMaxGNormRatioForWave`](/classes/wavevortextransform/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
    + [`u_max`](/classes/wavevortextransform/u_max.html) 
    + [`variables`](/classes/wavevortextransform/variables.html) Primary method for accessing the dynamical variables on the
    + [`variablesAtPosition`](/classes/wavevortextransform/variablesatposition.html) Primary method for accessing the dynamical variables on the
    + [`velocityField`](/classes/wavevortextransform/velocityfield.html) Return the velocity field, which is the sum of the gridded
    + [`version`](/classes/wavevortextransform/version.html) 
    + [`waveVortexCoefficientsAtTimeT`](/classes/wavevortextransform/wavevortexcoefficientsattimet.html) 
+ State Variables
  + Other
    + [`Fm`](/classes/wavevortextransform/fm.html) non-linear flux into Am
    + [`Fp`](/classes/wavevortextransform/fp.html) non-linear flux into Ap
    + [`eta`](/classes/wavevortextransform/eta.html) isopycnal deviation
    + [`p`](/classes/wavevortextransform/p.html) pressure anomaly
    + [`qgpv`](/classes/wavevortextransform/qgpv.html) quasigeostrophic potential vorticity
    + [`u`](/classes/wavevortextransform/u.html) x-component of the fluid velocity
    + [`v`](/classes/wavevortextransform/v.html) y-component of the fluid velocity
    + [`w`](/classes/wavevortextransform/w.html) z-component of the fluid velocity
+ External (non-gridded) modes
  + Other
    + [`addExternalWavesWithFrequencies`](/classes/wavevortextransform/addexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given wavenumber
    + [`addExternalWavesWithWavenumbers`](/classes/wavevortextransform/addexternalwaveswithwavenumbers.html) add external (non-gridded) waves with a given wavenumber
    + [`externalVariableFieldsAtTime`](/classes/wavevortextransform/externalvariablefieldsattime.html) Returns the external wave modes at the grid points.
    + [`externalVariablesAtTimePosition`](/classes/wavevortextransform/externalvariablesattimeposition.html) Returns the external wave modes at the grid points.
    + [`fillOutWaveSpectrum`](/classes/wavevortextransform/filloutwavespectrum.html) Add external waves to the model to fill out the spectrum
    + [`removeAllExternalWaves`](/classes/wavevortextransform/removeallexternalwaves.html) remove all external (non-gridded) waves
    + [`setExternalWavesWithFrequencies`](/classes/wavevortextransform/setexternalwaveswithfrequencies.html) set external (non-gridded) waves with a given frequency
    + [`setExternalWavesWithWavenumbers`](/classes/wavevortextransform/setexternalwaveswithwavenumbers.html) set external (non-gridded) waves with a given wavenumber
+ Internal
  + Other
    + [`addToVariableCache`](/classes/wavevortextransform/addtovariablecache.html) add variable to internal cache, in case it is needed again
    + [`clearVariableCache`](/classes/wavevortextransform/clearvariablecache.html) clear the internal cache
    + [`clearVariableCacheOfTimeDependentVariables`](/classes/wavevortextransform/clearvariablecacheoftimedependentvariables.html) clear the internal cache of variables that claim to be time dependent
    + [`defaultTransformDimensions`](/classes/wavevortextransform/defaulttransformdimensions.html) return array of TransformDimensions initialized by default
    + [`defaultTransformMethods`](/classes/wavevortextransform/defaulttransformmethods.html) return array of TransformAnnotations to annotate the methods
    + [`defaultTransformOperations`](/classes/wavevortextransform/defaulttransformoperations.html) return array of TransformOperation instances initialized by default
    + [`defaultTransformProperties`](/classes/wavevortextransform/defaulttransformproperties.html) return array of TransformProperty initialized by default
    + [`fetchFromVariableCache`](/classes/wavevortextransform/fetchfromvariablecache.html) retrieve a set of variables from the internal cache
    + [`performTransformOperation`](/classes/wavevortextransform/performtransformoperation.html) computes (runs) the operation
    + [`stateVariables`](/classes/wavevortextransform/statevariables.html) retrieve variables either from cache or by computation
+ Utility function
  + Metadata
    + [`addTransformDimension`](/classes/wavevortextransform/addtransformdimension.html) add a TransformDimension
    + [`addTransformProperty`](/classes/wavevortextransform/addtransformproperty.html) add a addTransformProperty
  + Other
    + [`checkHermitian`](/classes/wavevortextransform/checkhermitian.html) Check if the matrix is Hermitian. Report errors.
    + [`extractNonzeroWaveProperties`](/classes/wavevortextransform/extractnonzerowaveproperties.html) Takes a Hermitian matrix and returns the amplitude and phase of nonzero components
    + [`generateHermitianRandomMatrix`](/classes/wavevortextransform/generatehermitianrandommatrix.html) Generate a 3D matrix to be Hermitian, except at k=l=0
    + [`makeHermitian`](/classes/wavevortextransform/makehermitian.html) Forces a 3D matrix to be Hermitian
    + [`nyquistWavenumbers`](/classes/wavevortextransform/nyquistwavenumbers.html) Returns a matrix with 1s at the Nyquist frequencies.
    + [`redundantHermitianCoefficients`](/classes/wavevortextransform/redundanthermitiancoefficients.html) Returns a matrix with 1s at the 'redundant' hermiation indices.
+ Operations
  + Create new operations and variables
    + [`addTransformOperation`](/classes/wavevortextransform/addtransformoperation.html) add a addTransformProperty
  + differentiation
    + [`diffX`](/classes/wavevortextransform/diffx.html) differentiate a spatial variable in the x-direction
    + [`diffY`](/classes/wavevortextransform/diffy.html) differentiate a spatial variable in the y-direction
  + Transformations
    + [`transformFromSpatialDomainWithF`](/classes/wavevortextransform/transformfromspatialdomainwithf.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the F-modes
    + [`transformFromSpatialDomainWithG`](/classes/wavevortextransform/transformfromspatialdomainwithg.html) transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the G-modes
    + [`transformToSpatialDomainWithF`](/classes/wavevortextransform/transformtospatialdomainwithf.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes
    + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithfallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.
    + [`transformToSpatialDomainWithG`](/classes/wavevortextransform/transformtospatialdomainwithg.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes
    + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithgallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.
+ Masks
  + Other
    + [`maskForAliasedModes`](/classes/wavevortextransform/maskforaliasedmodes.html) returns a mask with locations of modes that will alias with a quadratic multiplication.
    + [`masksForAllFlowConstituents`](/classes/wavevortextransform/masksforallflowconstituents.html) Returns six 'masks' (matrices with 1s or 0s) indicating where the six
    + [`masksForFlowConstituents`](/classes/wavevortextransform/masksforflowconstituents.html) Returns a sets of 'masks' indicating where different solution types live in the Ap, Am, A0 matrices.
+ Validation and internal unit testing
  + Other
    + [`validateTransformationMatrices`](/classes/wavevortextransform/validatetransformationmatrices.html) used to confirm if $$S$$ and $$S^{-1}$$ are inverses
+ Write to file
  + Other
    + [`writeToFile`](/classes/wavevortextransform/writetofile.html) Output the `WaveVortexTransform` to file.


---