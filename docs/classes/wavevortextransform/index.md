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
+ Other
  + Other
    + [`A0N`](/classes/wavevortextransform/a0n.html) matrix coefficient that multiplies $$\tilde{\eta}$$ to compute $$A_0$$.
    + [`A0U`](/classes/wavevortextransform/a0u.html) matrix coefficient that multiplies $$\tilde{u}$$ to compute $$A_0$$.
    + [`A0V`](/classes/wavevortextransform/a0v.html) matrix coefficient that multiplies $$\tilde{v}$$ to compute $$A_0$$.
    + [`A0_HKE_factor`](/classes/wavevortextransform/a0_hke_factor.html) 
    + [`A0_PE_factor`](/classes/wavevortextransform/a0_pe_factor.html) 
    + [`A0_TE_factor`](/classes/wavevortextransform/a0_te_factor.html) 
    + [`A0t`](/classes/wavevortextransform/a0t.html) geostrophic coefficients at time (t-t0)
    + [`AddExternalWavesWithFrequencies`](/classes/wavevortextransform/addexternalwaveswithfrequencies.html) 
    + [`AddExternalWavesWithWavenumbers`](/classes/wavevortextransform/addexternalwaveswithwavenumbers.html) 
    + [`Am`](/classes/wavevortextransform/am.html) negative wave coefficients at reference time t0
    + [`AmN`](/classes/wavevortextransform/amn.html) matrix coefficient that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.
    + [`AmU`](/classes/wavevortextransform/amu.html) matrix coefficient that multiplies $$\tilde{u}$$ to compute $$A_m$$.
    + [`AmV`](/classes/wavevortextransform/amv.html) matrix coefficient that multiplies $$\tilde{v}$$ to compute $$A_m$$.
    + [`Amt`](/classes/wavevortextransform/amt.html) negative wave coefficients at time (t-t0)
    + [`Ap`](/classes/wavevortextransform/ap.html) positive wave coefficients at reference time t0
    + [`ApN`](/classes/wavevortextransform/apn.html) matrix coefficient that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.
    + [`ApV`](/classes/wavevortextransform/apv.html) matrix coefficient that multiplies $$\tilde{v}$$ to compute $$A_p$$.
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
    + [`FillOutWaveSpectrum`](/classes/wavevortextransform/filloutwavespectrum.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`Fm`](/classes/wavevortextransform/fm.html) non-linear flux into Am
    + [`Fp`](/classes/wavevortextransform/fp.html) non-linear flux into Ap
    + [`InitializeWithRandomFlowState`](/classes/wavevortextransform/initializewithrandomflowstate.html) 
    + [`Kh`](/classes/wavevortextransform/kh.html) 
    + [`MaskForAliasedModes`](/classes/wavevortextransform/maskforaliasedmodes.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`MasksForAllFlowConstituents`](/classes/wavevortextransform/masksforallflowconstituents.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`MasksForFlowConstituents`](/classes/wavevortextransform/masksforflowconstituents.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`MasksForFlowContinuents`](/classes/wavevortextransform/masksforflowcontinuents.html) 
    + [`N0`](/classes/wavevortextransform/n0.html) interior buoyancy frequency at the surface (z=0)
    + [`N2`](/classes/wavevortextransform/n2.html) buoyancy frequency of the mean density
    + [`NA0`](/classes/wavevortextransform/na0.html) matrix coefficient that multiplies $$\A_0$$ to compute $$\tilde{\eta}$$.
    + [`NAm`](/classes/wavevortextransform/nam.html) matrix coefficient that multiplies $$\A_m$$ to compute $$\tilde{\eta}$$.
    + [`NAp`](/classes/wavevortextransform/nap.html) matrix coefficient that multiplies $$\A_p$$ to compute $$\tilde{\eta}$$.
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
    + [`UA0`](/classes/wavevortextransform/ua0.html) matrix coefficient that multiplies $$\A_0$$ to compute $$\tilde{u}$$.
    + [`UAm`](/classes/wavevortextransform/uam.html) matrix coefficient that multiplies $$\A_m$$ to compute $$\tilde{u}$$.
    + [`UAp`](/classes/wavevortextransform/uap.html) matrix coefficient that multiplies $$\A_p$$ to compute $$\tilde{u}$$.
    + [`VA0`](/classes/wavevortextransform/va0.html) matrix coefficient that multiplies $$\A_0$$ to compute $$\tilde{v}$$.
    + [`VAm`](/classes/wavevortextransform/vam.html) matrix coefficient that multiplies $$\A_m$$ to compute $$\tilde{v}$$.
    + [`VAp`](/classes/wavevortextransform/vap.html) matrix coefficient that multiplies $$\A_p$$ to compute $$\tilde{v}$$.
    + [`ValidateTransformationMatrices`](/classes/wavevortextransform/validatetransformationmatrices.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`WAm`](/classes/wavevortextransform/wam.html) matrix coefficient that multiplies $$\A_m$$ to compute $$\tilde{w}$$.
    + [`WAp`](/classes/wavevortextransform/wap.html) matrix coefficient that multiplies $$\A_p$$ to compute $$\tilde{w}$$.
    + [`WaveVortexTransform`](/classes/wavevortextransform/wavevortextransform.html) These first properties are directly set on initialization
    + [`X`](/classes/wavevortextransform/x.html) 
    + [`Y`](/classes/wavevortextransform/y.html) 
    + [`Z`](/classes/wavevortextransform/z.html) 
    + [`addToVariableCache`](/classes/wavevortextransform/addtovariablecache.html) 
    + [`addTransformDimension`](/classes/wavevortextransform/addtransformdimension.html) 
    + [`addTransformOperation`](/classes/wavevortextransform/addtransformoperation.html) 
    + [`addTransformProperty`](/classes/wavevortextransform/addtransformproperty.html) 
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
    + [`fetchFromVariableCache`](/classes/wavevortextransform/fetchfromvariablecache.html) 
    + [`g`](/classes/wavevortextransform/g.html) 
    + [`generateHermitianRandomMatrix`](/classes/wavevortextransform/generatehermitianrandommatrix.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`generateRandomFlowState`](/classes/wavevortextransform/generaterandomflowstate.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`geostrophicEnergy`](/classes/wavevortextransform/geostrophicenergy.html) 
    + [`geostrophicEnergyBaroclinic`](/classes/wavevortextransform/geostrophicenergybaroclinic.html) total energy, geostrophic, baroclinic
    + [`geostrophicEnergyBarotropic`](/classes/wavevortextransform/geostrophicenergybarotropic.html) total energy, geostrophic, barotropic
    + [`grid`](/classes/wavevortextransform/grid.html) 
    + [`h`](/classes/wavevortextransform/h.html) equivalent depth of each mode
    + [`iOmega`](/classes/wavevortextransform/iomega.html) 
    + [`inertialEnergy`](/classes/wavevortextransform/inertialenergy.html) 
    + [`inertialEnergyBaroclinic`](/classes/wavevortextransform/inertialenergybaroclinic.html) total energy, inertial oscillations, baroclinic
    + [`inertialEnergyBarotropic`](/classes/wavevortextransform/inertialenergybarotropic.html) total energy, inertial oscillations, barotropic
    + [`initWithGMSpectrum`](/classes/wavevortextransform/initwithgmspectrum.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    + [`initWithSpectralFunction`](/classes/wavevortextransform/initwithspectralfunction.html) The GM2D_int function is used to assign variance to a given
    + [`internalWaveEnergyMinus`](/classes/wavevortextransform/internalwaveenergyminus.html) total energy, internal waves, minus
    + [`internalWaveEnergyPlus`](/classes/wavevortextransform/internalwaveenergyplus.html) total energy, internal waves, positive
    + [`interpolatedFieldAtPositionBadBoundaries`](/classes/wavevortextransform/interpolatedfieldatpositionbadboundaries.html) 
    + [`latitude`](/classes/wavevortextransform/latitude.html) latitude of the simulation
    + [`makeHermitian`](/classes/wavevortextransform/makehermitian.html) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    + [`transformFromFile`](/classes/wavevortextransform/transformfromfile.html) 
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
+ Transformation matrix coefficients
  + Other
    + [`ApU`](/classes/wavevortextransform/apu.html) matrix coefficient that multiplies $$\tilde{u}$$ to compute $$A_p$$.
+ Write to file
  + Other
    + [`writeToFile`](/classes/wavevortextransform/writetofile.html) Output the `WaveVortexTransform` to file.


---