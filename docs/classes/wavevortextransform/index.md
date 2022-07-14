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
  state of the ocean at a given instant in time (e.g., u, v, w, and
  rho). What makes the WaveVortexTransform subclasses special is that
  the state of the ocean is represented as energetically independent
  waves and geostrophic motions (vortices). These classes can be
  queried for other information, e.g., Ertel PV, relative vorticity,
  etc. Then we have $$\exp(i\omega t)$$
 
                        


## Topics
+ Initialization
  + Other
    + [`transformFromFile`](/classes/wavevortextransform/transformfromfile.html) Initialize a WaveVortexTransform instance from an existing file
+ Domain attributes
  + Grid
    + [`Lx`](/classes/wavevortextransform/lx.html) domain size in the x-direction
    + [`Ly`](/classes/wavevortextransform/ly.html) domain size in the y-direction
    + [`Lz`](/classes/wavevortextransform/lz.html) domain size in the z-direction
    + [`j`](/classes/wavevortextransform/j.html) vertical mode number
    + [`k`](/classes/wavevortextransform/k.html) wavenumber-coordinate dimension in the x-direction
    + [`l`](/classes/wavevortextransform/l.html) wavenumber-coordinate dimension in the y-direction
    + [`x`](/classes/wavevortextransform/x.html) x-coordinate dimension
    + [`y`](/classes/wavevortextransform/y.html) y-coordinate dimension
    + [`z`](/classes/wavevortextransform/z.html) z-coordinate dimension
  + Other
    + [`F0`](/classes/wavevortextransform/f0.html) non-linear flux into A0
    + [`f0`](/classes/wavevortextransform/f0.html) Coriolis parameter
    + [`inertialPeriod`](/classes/wavevortextransform/inertialperiod.html) inertial period
    + [`isBarotropic`](/classes/wavevortextransform/isbarotropic.html) Boolean indicating whether there is a single (equivalent barotropic) mode
    + [`kRadial`](/classes/wavevortextransform/kradial.html) isotropic wavenumber dimension
    + [`t`](/classes/wavevortextransform/t.html) time of observations
    + [`t0`](/classes/wavevortextransform/t0.html) reference time of Ap, Am, A0
+ Wave-vortex coefficients
  + Other
    + [`A0`](/classes/wavevortextransform/a0.html) geostrophic coefficients at reference time t0
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
  + Geostrophic Constituents
  + Inertia-Gravity Wave Constituents
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
    + [`A0_HKE_factor`](/classes/wavevortextransform/a0_hke_factor.html) 
    + [`A0_PE_factor`](/classes/wavevortextransform/a0_pe_factor.html) 
    + [`A0_TE_factor`](/classes/wavevortextransform/a0_te_factor.html) 
    + [`A0t`](/classes/wavevortextransform/a0t.html) geostrophic coefficients at time (t-t0)
    + [`Am`](/classes/wavevortextransform/am.html) negative wave coefficients at reference time t0
    + [`Amt`](/classes/wavevortextransform/amt.html) negative wave coefficients at time (t-t0)
    + [`Ap`](/classes/wavevortextransform/ap.html) positive wave coefficients at reference time t0
    + [`Apm_TE_factor`](/classes/wavevortextransform/apm_te_factor.html) These convert the coefficients to their depth integrated energies
    + [`Apt`](/classes/wavevortextransform/apt.html) positive wave coefficients at time (t-t0)
    + [`EnergeticsByWavenumberAndMode`](/classes/wavevortextransform/energeticsbywavenumberandmode.html) 
    + [`EnergyFluxAtTimeInitial`](/classes/wavevortextransform/energyfluxattimeinitial.html) 
    + [`EnergyFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/energyfluxforflowconstituentsattime.html) 
    + [`EnergyFluxForFlowConstituentsAtTimeInitial`](/classes/wavevortextransform/energyfluxforflowconstituentsattimeinitial.html) 
    + [`ExponentialFilter`](/classes/wavevortextransform/exponentialfilter.html) 
    + [`Fm`](/classes/wavevortextransform/fm.html) non-linear flux into Am
    + [`Fp`](/classes/wavevortextransform/fp.html) non-linear flux into Ap
    + [`InitializeWithRandomFlowState`](/classes/wavevortextransform/initializewithrandomflowstate.html) 
    + [`Kh`](/classes/wavevortextransform/kh.html) 
    + [`N0`](/classes/wavevortextransform/n0.html) interior buoyancy frequency at the surface (z=0)
    + [`N2`](/classes/wavevortextransform/n2.html) buoyancy frequency of the mean density
    + [`Nj`](/classes/wavevortextransform/nj.html) 
    + [`Nk`](/classes/wavevortextransform/nk.html) 
    + [`Nl`](/classes/wavevortextransform/nl.html) 
    + [`Nmax`](/classes/wavevortextransform/nmax.html) maximum buoyancy frequency
    + [`NonlinearFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/nonlinearfluxforflowconstituentsattime.html) Apply operator T_\omega---defined in (C2) in the manuscript
    + [`Nx`](/classes/wavevortextransform/nx.html) 
    + [`Ny`](/classes/wavevortextransform/ny.html) 
    + [`Nz`](/classes/wavevortextransform/nz.html) 
    + [`Omega`](/classes/wavevortextransform/omega.html) 
    + [`WaveVortexTransform`](/classes/wavevortextransform/wavevortextransform.html) These first properties are directly set on initialization
    + [`X`](/classes/wavevortextransform/x.html) 
    + [`Y`](/classes/wavevortextransform/y.html) 
    + [`Z`](/classes/wavevortextransform/z.html) 
    + [`addToVariableCache`](/classes/wavevortextransform/addtovariablecache.html) 
    + [`addTransformDimension`](/classes/wavevortextransform/addtransformdimension.html) 
    + [`addTransformOperation`](/classes/wavevortextransform/addtransformoperation.html) 
    + [`addTransformProperty`](/classes/wavevortextransform/addtransformproperty.html) 
    + [`buildTransformationMatrices`](/classes/wavevortextransform/buildtransformationmatrices.html) Build wavenumbers
    + [`clearVariableCache`](/classes/wavevortextransform/clearvariablecache.html) 
    + [`clearVariableCacheOfTimeDependentVariables`](/classes/wavevortextransform/clearvariablecacheoftimedependentvariables.html) 
    + [`dLnN2`](/classes/wavevortextransform/dlnn2.html) d/dz ln N2
    + [`defaultTransformDimensions`](/classes/wavevortextransform/defaulttransformdimensions.html) 
    + [`defaultTransformOperations`](/classes/wavevortextransform/defaulttransformoperations.html) 
    + [`defaultTransformProperties`](/classes/wavevortextransform/defaulttransformproperties.html) 
    + [`diffX`](/classes/wavevortextransform/diffx.html) 
    + [`diffY`](/classes/wavevortextransform/diffy.html) 
    + [`diffZF`](/classes/wavevortextransform/diffzf.html) 
    + [`diffZG`](/classes/wavevortextransform/diffzg.html) 
    + [`energyFlux`](/classes/wavevortextransform/energyflux.html) 
    + [`eta`](/classes/wavevortextransform/eta.html) isopycnal deviation
    + [`fetchFromVariableCache`](/classes/wavevortextransform/fetchfromvariablecache.html) 
    + [`g`](/classes/wavevortextransform/g.html) 
    + [`generateRandomFlowState`](/classes/wavevortextransform/generaterandomflowstate.html) Random flow state, separated out by solution type.
    + [`geostrophicEnergy`](/classes/wavevortextransform/geostrophicenergy.html) 
    + [`geostrophicEnergyBaroclinic`](/classes/wavevortextransform/geostrophicenergybaroclinic.html) total energy, geostrophic, baroclinic
    + [`geostrophicEnergyBarotropic`](/classes/wavevortextransform/geostrophicenergybarotropic.html) total energy, geostrophic, barotropic
    + [`grid`](/classes/wavevortextransform/grid.html) 
    + [`h`](/classes/wavevortextransform/h.html) equivalent depth of each mode
    + [`iOmega`](/classes/wavevortextransform/iomega.html) 
    + [`inertialEnergy`](/classes/wavevortextransform/inertialenergy.html) 
    + [`inertialEnergyBaroclinic`](/classes/wavevortextransform/inertialenergybaroclinic.html) total energy, inertial oscillations, baroclinic
    + [`inertialEnergyBarotropic`](/classes/wavevortextransform/inertialenergybarotropic.html) total energy, inertial oscillations, barotropic
    + [`internalWaveEnergyMinus`](/classes/wavevortextransform/internalwaveenergyminus.html) total energy, internal waves, minus
    + [`internalWaveEnergyPlus`](/classes/wavevortextransform/internalwaveenergyplus.html) total energy, internal waves, positive
    + [`interpolatedFieldAtPositionBadBoundaries`](/classes/wavevortextransform/interpolatedfieldatpositionbadboundaries.html) 
    + [`latitude`](/classes/wavevortextransform/latitude.html) latitude of the simulation
    + [`nonlinearFlux`](/classes/wavevortextransform/nonlinearflux.html) 
    + [`offgridModes`](/classes/wavevortextransform/offgridmodes.html) offgridModes -  subclass should initialize
    + [`ongridModes`](/classes/wavevortextransform/ongridmodes.html) ongridModes -  This is a cached copy
    + [`p`](/classes/wavevortextransform/p.html) pressure anomaly
    + [`performTransformOperation`](/classes/wavevortextransform/performtransformoperation.html) 
    + [`qgpv`](/classes/wavevortextransform/qgpv.html) quasigeostrophic potential vorticity
    + [`radialWavenumberAxis`](/classes/wavevortextransform/radialwavenumberaxis.html) Create a reasonable wavenumber axis
    + [`rebuildTransformationMatrices`](/classes/wavevortextransform/rebuildtransformationmatrices.html) 
    + [`rho0`](/classes/wavevortextransform/rho0.html) mean density at the surface (z=0)
    + [`rhobar`](/classes/wavevortextransform/rhobar.html) mean density
    + [`spectralVanishingViscosityFilter`](/classes/wavevortextransform/spectralvanishingviscosityfilter.html) Builds the spectral vanishing viscosity operator
    + [`stateVariableWithName`](/classes/wavevortextransform/statevariablewithname.html) 
    + [`stateVariables`](/classes/wavevortextransform/statevariables.html) 
    + [`summarizeEnergyContent`](/classes/wavevortextransform/summarizeenergycontent.html) 
    + [`timeDependentStateVariables`](/classes/wavevortextransform/timedependentstatevariables.html) 
    + [`totalEnergy`](/classes/wavevortextransform/totalenergy.html) energy = self.inertialEnergy + self.waveEnergy + self.geostrophicEnergy;
    + [`totalEnergySpatiallyIntegrated`](/classes/wavevortextransform/totalenergyspatiallyintegrated.html) 
    + [`totalHydrostaticEnergy`](/classes/wavevortextransform/totalhydrostaticenergy.html) 
    + [`transformDimensionWithName`](/classes/wavevortextransform/transformdimensionwithname.html) 
    + [`transformFromSpatialDomainWithF`](/classes/wavevortextransform/transformfromspatialdomainwithf.html) 
    + [`transformFromSpatialDomainWithG`](/classes/wavevortextransform/transformfromspatialdomainwithg.html) 
    + [`transformOperationWithName`](/classes/wavevortextransform/transformoperationwithname.html) 
    + [`transformPropertyWithName`](/classes/wavevortextransform/transformpropertywithname.html) 
    + [`transformToRadialWavenumber`](/classes/wavevortextransform/transformtoradialwavenumber.html) 
    + [`transformToSpatialDomainWithF`](/classes/wavevortextransform/transformtospatialdomainwithf.html) 
    + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithfallderivatives.html) 
    + [`transformToSpatialDomainWithG`](/classes/wavevortextransform/transformtospatialdomainwithg.html) 
    + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithgallderivatives.html) 
    + [`transformUVEtaToWaveVortex`](/classes/wavevortextransform/transformuvetatowavevortex.html) This is the 'S^{-1}' operator (C5) in the manuscript
    + [`transformWaveVortexToUVWEta`](/classes/wavevortextransform/transformwavevortextouvweta.html) 
    + [`transformWithDoubleResolution`](/classes/wavevortextransform/transformwithdoubleresolution.html) 
    + [`transformWithResolution`](/classes/wavevortextransform/transformwithresolution.html) 
    + [`u`](/classes/wavevortextransform/u.html) x-component of the fluid velocity
    + [`uMaxGNormRatioForWave`](/classes/wavevortextransform/umaxgnormratioforwave.html) Needed to add and remove internal waves from the model
    + [`u_max`](/classes/wavevortextransform/u_max.html) 
    + [`v`](/classes/wavevortextransform/v.html) y-component of the fluid velocity
    + [`variables`](/classes/wavevortextransform/variables.html) Primary method for accessing the dynamical variables on the
    + [`variablesAtPosition`](/classes/wavevortextransform/variablesatposition.html) Primary method for accessing the dynamical variables on the
    + [`velocityField`](/classes/wavevortextransform/velocityfield.html) Return the velocity field, which is the sum of the gridded
    + [`version`](/classes/wavevortextransform/version.html) 
    + [`w`](/classes/wavevortextransform/w.html) z-component of the fluid velocity
    + [`waveEnergy`](/classes/wavevortextransform/waveenergy.html) 
    + [`waveVortexCoefficientsAtTimeT`](/classes/wavevortextransform/wavevortexcoefficientsattimet.html) 
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
+ Utility function
  + Other
    + [`checkHermitian`](/classes/wavevortextransform/checkhermitian.html) Check if the matrix is Hermitian. Report errors.
    + [`extractNonzeroWaveProperties`](/classes/wavevortextransform/extractnonzerowaveproperties.html) Takes a Hermitian matrix and returns the amplitude and phase of nonzero components
    + [`generateHermitianRandomMatrix`](/classes/wavevortextransform/generatehermitianrandommatrix.html) Generate a 3D matrix to be Hermitian, except at k=l=0
    + [`makeHermitian`](/classes/wavevortextransform/makehermitian.html) Forces a 3D matrix to be Hermitian
    + [`nyquistWavenumbers`](/classes/wavevortextransform/nyquistwavenumbers.html) Returns a matrix with 1s at the Nyquist frequencies.
    + [`redundantHermitianCoefficients`](/classes/wavevortextransform/redundanthermitiancoefficients.html) Returns a matrix with 1s at the 'redundant' hermiation indices.
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