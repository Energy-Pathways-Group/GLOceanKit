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
+ Domain attributes
  + Other
    + [`f0`](/classes/wavevortextransform/f0.html) geostrophic coefficients at reference time t0
    + [`inertialPeriod`](/classes/wavevortextransform/inertialperiod.html) 
    + [`isBarotropic`](/classes/wavevortextransform/isbarotropic.html) 
    + [`j`](/classes/wavevortextransform/j.html) 
    + [`k`](/classes/wavevortextransform/k.html) 
    + [`l`](/classes/wavevortextransform/l.html) 
    + [`x`](/classes/wavevortextransform/x.html) 
    + [`y`](/classes/wavevortextransform/y.html) geostrophic coefficients at time (t-t0)
+ Other
  + Other
    + [`A0`](/classes/wavevortextransform/a0.html) geostrophic coefficients at reference time t0
    + [`A0N`](/classes/wavevortextransform/a0n.html) 
    + [`A0U`](/classes/wavevortextransform/a0u.html) 
    + [`A0V`](/classes/wavevortextransform/a0v.html) 
    + [`A0_HKE_factor`](/classes/wavevortextransform/a0_hke_factor.html) 
    + [`A0_PE_factor`](/classes/wavevortextransform/a0_pe_factor.html) 
    + [`A0_TE_factor`](/classes/wavevortextransform/a0_te_factor.html) 
    + [`A0t`](/classes/wavevortextransform/a0t.html) geostrophic coefficients at time (t-t0)
    + [`AddExternalWavesWithFrequencies`](/classes/wavevortextransform/addexternalwaveswithfrequencies.html) 
    + [`AddExternalWavesWithWavenumbers`](/classes/wavevortextransform/addexternalwaveswithwavenumbers.html) 
    + [`Am`](/classes/wavevortextransform/am.html) negative wave coefficients at reference time t0
    + [`AmN`](/classes/wavevortextransform/amn.html) 
    + [`AmU`](/classes/wavevortextransform/amu.html) 
    + [`AmV`](/classes/wavevortextransform/amv.html) 
    + [`Amt`](/classes/wavevortextransform/amt.html) negative wave coefficients at time (t-t0)
    + [`Ap`](/classes/wavevortextransform/ap.html) positive wave coefficients at reference time t0
    + [`ApN`](/classes/wavevortextransform/apn.html) 
    + [`ApU`](/classes/wavevortextransform/apu.html) 
    + [`ApV`](/classes/wavevortextransform/apv.html) 
    + [`Apm_TE_factor`](/classes/wavevortextransform/apm_te_factor.html) These convert the coefficients to their depth integrated energies
    + [`Apt`](/classes/wavevortextransform/apt.html) positive wave coefficients at time (t-t0)
    + [`EnergeticsByWavenumberAndMode`](/classes/wavevortextransform/energeticsbywavenumberandmode.html) 
    + [`EnergyFluxAtTimeInitial`](/classes/wavevortextransform/energyfluxattimeinitial.html) 
    + [`EnergyFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/energyfluxforflowconstituentsattime.html) 
    + [`EnergyFluxForFlowConstituentsAtTimeInitial`](/classes/wavevortextransform/energyfluxforflowconstituentsattimeinitial.html) 
    + [`ExponentialFilter`](/classes/wavevortextransform/exponentialfilter.html) 
    + [`ExternalVariableFieldsAtTime`](/classes/wavevortextransform/externalvariablefieldsattime.html) Returns the external wave modes at the grid points.
    + [`ExternalVariablesAtTimePosition`](/classes/wavevortextransform/externalvariablesattimeposition.html) 
    + [`ExtractNonzeroWaveProperties`](/classes/wavevortextransform/extractnonzerowaveproperties.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`F0`](/classes/wavevortextransform/f0.html) non-linear flux into A0
    + [`FillOutWaveSpectrum`](/classes/wavevortextransform/filloutwavespectrum.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`Fm`](/classes/wavevortextransform/fm.html) non-linear flux into Am
    + [`Fp`](/classes/wavevortextransform/fp.html) non-linear flux into Ap
    + [`InitializeWithRandomFlowState`](/classes/wavevortextransform/initializewithrandomflowstate.html) 
    + [`Kh`](/classes/wavevortextransform/kh.html) 
    + [`Lx`](/classes/wavevortextransform/lx.html) domain size in the x-direction
    + [`Ly`](/classes/wavevortextransform/ly.html) domain size in the y-direction
    + [`Lz`](/classes/wavevortextransform/lz.html) domain size in the z-direction
    + [`MaskForAliasedModes`](/classes/wavevortextransform/maskforaliasedmodes.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`MasksForAllFlowConstituents`](/classes/wavevortextransform/masksforallflowconstituents.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`MasksForFlowConstituents`](/classes/wavevortextransform/masksforflowconstituents.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`MasksForFlowContinuents`](/classes/wavevortextransform/masksforflowcontinuents.html) 
    + [`N0`](/classes/wavevortextransform/n0.html) interior buoyancy frequency at the surface (z=0)
    + [`N2`](/classes/wavevortextransform/n2.html) buoyancy frequency of the mean density
    + [`NA0`](/classes/wavevortextransform/na0.html) 
    + [`NAm`](/classes/wavevortextransform/nam.html) 
    + [`NAp`](/classes/wavevortextransform/nap.html) 
    + [`Nj`](/classes/wavevortextransform/nj.html) 
    + [`Nk`](/classes/wavevortextransform/nk.html) 
    + [`Nl`](/classes/wavevortextransform/nl.html) 
    + [`Nmax`](/classes/wavevortextransform/nmax.html) maximum buoyancy frequency
    + [`NonlinearFluxForFlowConstituentsAtTime`](/classes/wavevortextransform/nonlinearfluxforflowconstituentsattime.html) Apply operator T_\omega---defined in (C2) in the manuscript
    + [`Nx`](/classes/wavevortextransform/nx.html) 
    + [`Ny`](/classes/wavevortextransform/ny.html) 
    + [`NyquistWavenumbers`](/classes/wavevortextransform/nyquistwavenumbers.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`Nz`](/classes/wavevortextransform/nz.html) 
    + [`Omega`](/classes/wavevortextransform/omega.html) 
    + [`RedundantHermitianCoefficients`](/classes/wavevortextransform/redundanthermitiancoefficients.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`RemoveAllExternalWaves`](/classes/wavevortextransform/removeallexternalwaves.html) 
    + [`SetExternalWavesWithFrequencies`](/classes/wavevortextransform/setexternalwaveswithfrequencies.html) 
    + [`SetExternalWavesWithWavenumbers`](/classes/wavevortextransform/setexternalwaveswithwavenumbers.html) 
    + [`UA0`](/classes/wavevortextransform/ua0.html) 
    + [`UAm`](/classes/wavevortextransform/uam.html) 
    + [`UAp`](/classes/wavevortextransform/uap.html) 
    + [`VA0`](/classes/wavevortextransform/va0.html) 
    + [`VAm`](/classes/wavevortextransform/vam.html) 
    + [`VAp`](/classes/wavevortextransform/vap.html) 
    + [`ValidateTransformationMatrices`](/classes/wavevortextransform/validatetransformationmatrices.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`WAm`](/classes/wavevortextransform/wam.html) 
    + [`WAp`](/classes/wavevortextransform/wap.html) 
    + [`WaveVortexTransform`](/classes/wavevortextransform/wavevortextransform.html) These first properties are directly set on initialization
    + [`X`](/classes/wavevortextransform/x.html) 
    + [`Y`](/classes/wavevortextransform/y.html) 
    + [`Z`](/classes/wavevortextransform/z.html) 
    + [`addGeostrophicStreamfunction`](/classes/wavevortextransform/addgeostrophicstreamfunction.html) 
    + [`addInertialMotions`](/classes/wavevortextransform/addinertialmotions.html) 
    + [`addToVariableCache`](/classes/wavevortextransform/addtovariablecache.html) 
    + [`addTransformDimension`](/classes/wavevortextransform/addtransformdimension.html) 
    + [`addTransformOperation`](/classes/wavevortextransform/addtransformoperation.html) 
    + [`addTransformProperty`](/classes/wavevortextransform/addtransformproperty.html) 
    + [`addWaveModes`](/classes/wavevortextransform/addwavemodes.html) 
    + [`buildTransformationMatrices`](/classes/wavevortextransform/buildtransformationmatrices.html) Build wavenumbers
    + [`checkHermitian`](/classes/wavevortextransform/checkhermitian.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    + [`fetchFromVariableCache`](/classes/wavevortextransform/fetchfromvariablecache.html) Coriolis parameter (radians/s)
    + [`g`](/classes/wavevortextransform/g.html) 
    + [`generateHermitianRandomMatrix`](/classes/wavevortextransform/generatehermitianrandommatrix.html) 
    + [`generateRandomFlowState`](/classes/wavevortextransform/generaterandomflowstate.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`geostrophicEnergy`](/classes/wavevortextransform/geostrophicenergy.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`geostrophicEnergyBaroclinic`](/classes/wavevortextransform/geostrophicenergybaroclinic.html) 
    + [`geostrophicEnergyBarotropic`](/classes/wavevortextransform/geostrophicenergybarotropic.html) total energy, geostrophic, baroclinic
    + [`grid`](/classes/wavevortextransform/grid.html) total energy, geostrophic, barotropic
    + [`h`](/classes/wavevortextransform/h.html) 
    + [`iOmega`](/classes/wavevortextransform/iomega.html) equivalent depth of each mode
    + [`inertialEnergy`](/classes/wavevortextransform/inertialenergy.html) 
    + [`inertialEnergyBaroclinic`](/classes/wavevortextransform/inertialenergybaroclinic.html) 
    + [`inertialEnergyBarotropic`](/classes/wavevortextransform/inertialenergybarotropic.html) total energy, inertial oscillations, baroclinic
    + [`initWithGMSpectrum`](/classes/wavevortextransform/initwithgmspectrum.html) total energy, inertial oscillations, barotropic
    + [`initWithGeostrophicStreamfunction`](/classes/wavevortextransform/initwithgeostrophicstreamfunction.html) Inertial period (s)
    + [`initWithInertialMotions`](/classes/wavevortextransform/initwithinertialmotions.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`initWithSpectralFunction`](/classes/wavevortextransform/initwithspectralfunction.html) 
    + [`initWithWaveModes`](/classes/wavevortextransform/initwithwavemodes.html) 
    + [`internalWaveEnergyMinus`](/classes/wavevortextransform/internalwaveenergyminus.html) The GM2D_int function is used to assign variance to a given
    + [`internalWaveEnergyPlus`](/classes/wavevortextransform/internalwaveenergyplus.html) 
    + [`interpolatedFieldAtPositionBadBoundaries`](/classes/wavevortextransform/interpolatedfieldatpositionbadboundaries.html) total energy, internal waves, minus
    + [`kRadial`](/classes/wavevortextransform/kradial.html) total energy, internal waves, positive
    + [`latitude`](/classes/wavevortextransform/latitude.html) 
    + [`makeHermitian`](/classes/wavevortextransform/makehermitian.html) Boolean indicating whether there is a single (equivalent barotropic) mode
    + [`nonlinearFlux`](/classes/wavevortextransform/nonlinearflux.html) vertical mode number
    + [`offgridModes`](/classes/wavevortextransform/offgridmodes.html) wavenumber-coordinate dimension in the x-direction
    + [`ongridModes`](/classes/wavevortextransform/ongridmodes.html) isotropic wavenumber dimension
    + [`p`](/classes/wavevortextransform/p.html) wavenumber-coordinate dimension in the y-direction
    + [`performTransformOperation`](/classes/wavevortextransform/performtransformoperation.html) latitude of the simulation
    + [`qgpv`](/classes/wavevortextransform/qgpv.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`radialWavenumberAxis`](/classes/wavevortextransform/radialwavenumberaxis.html) 
    + [`rebuildTransformationMatrices`](/classes/wavevortextransform/rebuildtransformationmatrices.html) offgridModes -  subclass should initialize
    + [`removeAllGeostrophicMotions`](/classes/wavevortextransform/removeallgeostrophicmotions.html) ongridModes -  This is a cached copy
    + [`removeAllInertialMotions`](/classes/wavevortextransform/removeallinertialmotions.html) pressure anomaly
    + [`removeAllWaves`](/classes/wavevortextransform/removeallwaves.html) 
    + [`rho0`](/classes/wavevortextransform/rho0.html) quasigeostrophic potential vorticity
    + [`rhobar`](/classes/wavevortextransform/rhobar.html) Create a reasonable wavenumber axis
    + [`setGeostrophicStreamfunction`](/classes/wavevortextransform/setgeostrophicstreamfunction.html) 
    + [`setInertialMotions`](/classes/wavevortextransform/setinertialmotions.html) 
    + [`setWaveModes`](/classes/wavevortextransform/setwavemodes.html) 
    + [`spectralVanishingViscosityFilter`](/classes/wavevortextransform/spectralvanishingviscosityfilter.html) 
    + [`stateVariableWithName`](/classes/wavevortextransform/statevariablewithname.html) mean density at the surface (z=0)
    + [`stateVariables`](/classes/wavevortextransform/statevariables.html) mean density
    + [`summarizeEnergyContent`](/classes/wavevortextransform/summarizeenergycontent.html) 
    + [`t`](/classes/wavevortextransform/t.html) 
    + [`t0`](/classes/wavevortextransform/t0.html) 
    + [`timeDependentStateVariables`](/classes/wavevortextransform/timedependentstatevariables.html) Builds the spectral vanishing viscosity operator
    + [`totalEnergy`](/classes/wavevortextransform/totalenergy.html) 
    + [`totalEnergySpatiallyIntegrated`](/classes/wavevortextransform/totalenergyspatiallyintegrated.html) 
    + [`totalHydrostaticEnergy`](/classes/wavevortextransform/totalhydrostaticenergy.html) 
    + [`transformDimensionWithName`](/classes/wavevortextransform/transformdimensionwithname.html) time of observations
    + [`transformFromFile`](/classes/wavevortextransform/transformfromfile.html) reference time of Ap, Am, A0
    + [`transformFromSpatialDomainWithF`](/classes/wavevortextransform/transformfromspatialdomainwithf.html) 
    + [`transformFromSpatialDomainWithG`](/classes/wavevortextransform/transformfromspatialdomainwithg.html) energy = self.inertialEnergy + self.waveEnergy + self.geostrophicEnergy;
    + [`transformOperationWithName`](/classes/wavevortextransform/transformoperationwithname.html) 
    + [`transformPropertyWithName`](/classes/wavevortextransform/transformpropertywithname.html) 
    + [`transformToRadialWavenumber`](/classes/wavevortextransform/transformtoradialwavenumber.html) 
    + [`transformToSpatialDomainWithF`](/classes/wavevortextransform/transformtospatialdomainwithf.html) 
    + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithfallderivatives.html) 
    + [`transformToSpatialDomainWithG`](/classes/wavevortextransform/transformtospatialdomainwithg.html) 
    + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wavevortextransform/transformtospatialdomainwithgallderivatives.html) 
    + [`transformUVEtaToWaveVortex`](/classes/wavevortextransform/transformuvetatowavevortex.html) 
    + [`transformWaveVortexToUVWEta`](/classes/wavevortextransform/transformwavevortextouvweta.html) 
    + [`transformWithDoubleResolution`](/classes/wavevortextransform/transformwithdoubleresolution.html) 
    + [`transformWithResolution`](/classes/wavevortextransform/transformwithresolution.html) 
    + [`u`](/classes/wavevortextransform/u.html) 
    + [`uMaxGNormRatioForWave`](/classes/wavevortextransform/umaxgnormratioforwave.html) 
    + [`u_max`](/classes/wavevortextransform/u_max.html) This is the 'S^{-1}' operator (C5) in the manuscript
    + [`v`](/classes/wavevortextransform/v.html) 
    + [`variables`](/classes/wavevortextransform/variables.html) 
    + [`variablesAtPosition`](/classes/wavevortextransform/variablesatposition.html) 
    + [`velocityField`](/classes/wavevortextransform/velocityfield.html) x-component of the fluid velocity
    + [`version`](/classes/wavevortextransform/version.html) Needed to add and remove internal waves from the model
    + [`w`](/classes/wavevortextransform/w.html) 
    + [`waveCoefficientsFromWaveModes`](/classes/wavevortextransform/wavecoefficientsfromwavemodes.html) y-component of the fluid velocity
    + [`waveEnergy`](/classes/wavevortextransform/waveenergy.html) Primary method for accessing the dynamical variables on the
    + [`waveModesFromWaveCoefficients`](/classes/wavevortextransform/wavemodesfromwavecoefficients.html) Primary method for accessing the dynamical variables on the
    + [`waveVortexCoefficientsAtTimeT`](/classes/wavevortextransform/wavevortexcoefficientsattimet.html) Return the velocity field, which is the sum of the gridded
    + [`writeToFile`](/classes/wavevortextransform/writetofile.html) 
    + [`z`](/classes/wavevortextransform/z.html) z-component of the fluid velocity


---