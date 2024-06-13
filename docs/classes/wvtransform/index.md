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
+ Domain attributes
  + [`Omega`](/classes/wvtransform/omega.html) frequency of oscillation of the linear waves
  + [`f`](/classes/wvtransform/f.html) Coriolis parameter
  + [`g`](/classes/wvtransform/g.html) gravity of Earth
  + [`inertialPeriod`](/classes/wvtransform/inertialperiod.html) inertial period
  + [`isBarotropic`](/classes/wvtransform/isbarotropic.html) Boolean indicating whether there is a single (equivalent barotropic) mode
  + [`latitude`](/classes/wvtransform/latitude.html) central latitude of the simulation
  + [`rho0`](/classes/wvtransform/rho0.html) density of $$\rho_\textrm{nm}$$ at the surface (z=0)
  + [`t`](/classes/wvtransform/t.html) time coordinate
  + [`t0`](/classes/wvtransform/t0.html) reference time of Ap, Am, A0
  + Grid
    + [`shouldAntialias`](/classes/wvtransform/shouldantialias.html) whether antialiasing is enabled
    + Spectral
      + [`J`](/classes/wvtransform/j_.html) j-coordinate matrix
      + [`K`](/classes/wvtransform/k_.html) k-coordinate matrix
      + [`Kh`](/classes/wvtransform/kh.html) horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$
      + [`L`](/classes/wvtransform/l_.html) l-coordinate matrix
      + [`Nj`](/classes/wvtransform/nj.html) points in the j-coordinate, `length(z)`
      + [`Nkl`](/classes/wvtransform/nkl.html) points in the kl-coordinate, `length(k)`
      + [`j`](/classes/wvtransform/j.html) vertical mode number
      + [`k`](/classes/wvtransform/k.html) wavenumber coordinate in the x-direction
      + [`kRadial`](/classes/wvtransform/kradial.html) isotropic wavenumber dimension
      + [`kljGrid`](/classes/wvtransform/kljgrid.html) returns the K, L, J coordinate matrices
      + [`l`](/classes/wvtransform/l.html) wavenumber coordinate in the y-direction
      + [`spectralMatrixSize`](/classes/wvtransform/spectralmatrixsize.html) returns the size of any spectral matrix, e.g., Ap, Am, A0
    + Spatial
      + [`x`](/classes/wvtransform/x.html) x coordinate
      + [`y`](/classes/wvtransform/y.html) y coordinate
      + [`z`](/classes/wvtransform/z.html) z coordinate
      + [`Lx`](/classes/wvtransform/lx.html) domain size in the x-direction
      + [`Ly`](/classes/wvtransform/ly.html) domain size in the y-direction
      + [`Lz`](/classes/wvtransform/lz.html) domain size in the z-direction
      + [`Nx`](/classes/wvtransform/nx.html) points in the x-coordinate, `length(x)`
      + [`Ny`](/classes/wvtransform/ny.html) points in the y-coordinate, `length(y)`
      + [`Nz`](/classes/wvtransform/nz.html) points in the z-coordinate, `length(z)`
      + [`X`](/classes/wvtransform/x_.html) x-coordinate matrix
      + [`Y`](/classes/wvtransform/y_.html) y-coordinate matrix
      + [`Z`](/classes/wvtransform/z_.html) z-coordinate matrix
      + [`spatialMatrixSize`](/classes/wvtransform/spatialmatrixsize.html) returns the size of all real-valued field variables
      + [`xyzGrid`](/classes/wvtransform/xyzgrid.html) returns the X, Y, Z coordinate matrices
+ Wave-vortex coefficients
  + at time $$t$$
    + [`A0t`](/classes/wvtransform/a0t.html) geostrophic coefficients at time t
    + [`Amt`](/classes/wvtransform/amt.html) negative wave coefficients at time t
    + [`Apt`](/classes/wvtransform/apt.html) positive wave coefficients at time t
+ Initial Conditions
  + [`addRandomFlow`](/classes/wvtransform/addrandomflow.html) add randomized flow to the existing state
  + [`addUVEta`](/classes/wvtransform/adduveta.html) add $$(u,v,\eta)$$ to the existing values
  + [`initFromNetCDFFile`](/classes/wvtransform/initfromnetcdffile.html) initialize the flow from a NetCDF file
  + [`initWithRandomFlow`](/classes/wvtransform/initwithrandomflow.html) initialize with a random flow state
  + [`initWithUVEta`](/classes/wvtransform/initwithuveta.html) initialize with fluid variables $$(u,v,\eta)$$
  + [`initWithUVRho`](/classes/wvtransform/initwithuvrho.html) initialize with fluid variables $$(u,v,\rho)$$
  + Waves
    + [`removeAll`](/classes/wvtransform/removeall.html) removes all energy from the model
+ Energetics
  + [`summarizeEnergyContent`](/classes/wvtransform/summarizeenergycontent.html) displays a summary of the energy content of the fluid
  + [`summarizeModeEnergy`](/classes/wvtransform/summarizemodeenergy.html) List the most energetic modes
  + Multiplicative factors
    + [`A0_TE_factor`](/classes/wvtransform/a0_te_factor.html) multiplicative factor that multiplies $$A_0^2$$ to compute total energy.
    + [`Apm_TE_factor`](/classes/wvtransform/apm_te_factor.html) multiplicative factor that multiplies $$A_\pm^2$$ to compute total energy.
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
+ Potential Vorticity & Enstrophy
  + Multiplicative factors
    + [`A0_QGPV_factor`](/classes/wvtransform/a0_qgpv_factor.html) multiplicative factor that multiplies $$A_0$$ to compute quasigeostrophic potential vorticity (QGPV).
    + [`A0_TZ_factor`](/classes/wvtransform/a0_tz_factor.html) multiplicative factor that multiplies $$A_0^2$$ to compute quasigeostrophic enstrophy.
+ State Variables
  + [`F0`](/classes/wvtransform/f0.html) non-linear flux into A0
  + [`Fm`](/classes/wvtransform/fm.html) non-linear flux into Am
  + [`Fp`](/classes/wvtransform/fp.html) non-linear flux into Ap
  + [`rho_e`](/classes/wvtransform/rho_e.html) excess density
  + [`rho_total`](/classes/wvtransform/rho_total.html) total potential density
  + [`ssh`](/classes/wvtransform/ssh.html) sea-surface height
  + [`ssu`](/classes/wvtransform/ssu.html) x-component of the fluid velocity at the surface
  + [`ssv`](/classes/wvtransform/ssv.html) y-component of the fluid velocity at the surface
  + [`uvMax`](/classes/wvtransform/uvmax.html) max horizontal fluid speed
  + [`wMax`](/classes/wvtransform/wmax.html) max vertical fluid speed
  + [`zeta_z`](/classes/wvtransform/zeta_z.html) vertical component of relative vorticity
+ Internal
  + [`WVTransform`](/classes/wvtransform/wvtransform.html) initialize a WVTransform instance
  + [`addToVariableCache`](/classes/wvtransform/addtovariablecache.html) add variable to internal cache, in case it is needed again
  + [`clearVariableCache`](/classes/wvtransform/clearvariablecache.html) clear the internal cache
  + [`clearVariableCacheOfTimeDependentVariables`](/classes/wvtransform/clearvariablecacheoftimedependentvariables.html) clear the internal cache of variables that claim to be time dependent
  + [`fetchFromVariableCache`](/classes/wvtransform/fetchfromvariablecache.html) retrieve a set of variables from the internal cache
  + [`performOperation`](/classes/wvtransform/performoperation.html) computes (runs) the operation
  + [`performOperationWithName`](/classes/wvtransform/performoperationwithname.html) computes (runs) the operation
  + [`stateVariables`](/classes/wvtransform/statevariables.html) retrieve variables either from cache or by computation
+ Metadata
  + Dimensions
    + [`addDimensionAnnotations`](/classes/wvtransform/adddimensionannotations.html) add one or more WVDimensions
    + [`dimensionAnnotationWithName`](/classes/wvtransform/dimensionannotationwithname.html) retrieve a WVDimension by name
+ Flow components
  + [`addFlowComponent`](/classes/wvtransform/addflowcomponent.html) add a flow component
  + [`addPrimaryFlowComponent`](/classes/wvtransform/addprimaryflowcomponent.html) add a primary flow component, automatically added to the flow
  + [`flowComponent`](/classes/wvtransform/flowcomponent.html) retrieve a WVFlowComponent by name
  + [`primaryFlowComponent`](/classes/wvtransform/primaryflowcomponent.html) retrieve a WVPrimaryFlowComponent by name
+ Utility function
  + [`spectralVariableWithResolution`](/classes/wvtransform/spectralvariablewithresolution.html) create a new variable with different resolution
  + Metadata
    + [`addOperation`](/classes/wvtransform/addoperation.html) add a WVOperation
    + [`addPropertyAnnotations`](/classes/wvtransform/addpropertyannotations.html) add a property annotation
    + [`addVariableAnnotations`](/classes/wvtransform/addvariableannotations.html) add a variable annotation
    + [`operationWithName`](/classes/wvtransform/operationwithname.html) retrieve a WVOperation by name
    + [`propertyAnnotationWithName`](/classes/wvtransform/propertyannotationwithname.html) retrieve a WVPropertyAnnotation by name
    + [`removeOperation`](/classes/wvtransform/removeoperation.html) remove an existing WVOperation
    + [`removeVariableAnnotations`](/classes/wvtransform/removevariableannotations.html) add a variable annotation
    + [`variableAnnotationWithName`](/classes/wvtransform/variableannotationwithname.html) retrieve a WVVariableAnnotation by name
    + [`variableNames`](/classes/wvtransform/variablenames.html) retrieve the names of all available variables
+ Write to file
  + [`concatenateVariablesAlongTimeDimension`](/classes/wvtransform/concatenatevariablesalongtimedimension.html) Concatenate variables along the time dimension
  + [`createNetCDFFileForTimeStepOutput`](/classes/wvtransform/createnetcdffilefortimestepoutput.html) Output the `WVTransform` to file with variable time dimension
  + [`writeToFile`](/classes/wvtransform/writetofile.html) Output the `WVTransform` to file.
+ Operations
  + Transformations
    + [`convertFromWavenumberToFrequency`](/classes/wvtransform/convertfromwavenumbertofrequency.html) Summary
    + [`transformFromSpatialDomainWithFg`](/classes/wvtransform/transformfromspatialdomainwithfg.html) transforms from the spatial domain (z,:,:) to the spectral domain (j,:,:) using the geostrophic F-modes
    + [`transformFromSpatialDomainWithFio`](/classes/wvtransform/transformfromspatialdomainwithfio.html) transforms from the spatial domain (z,:,:) to the spectral domain (j,:,:) using the inertial oscillation F-modes
    + [`transformFromSpatialDomainWithGg`](/classes/wvtransform/transformfromspatialdomainwithgg.html) transforms from the spatial domain (z,:,:) to the spectral domain (j,:,:) using the geostrophic G-modes
    + [`transformToKLAxes`](/classes/wvtransform/transformtoklaxes.html) transforms in the spectral domain from (j,kl) to (kAxis,lAxis,j)
    + [`transformToRadialWavenumber`](/classes/wvtransform/transformtoradialwavenumber.html) transforms in the spectral domain from (j,kl) to (j,kRadial)
    + [`transformToSpatialDomainWithF`](/classes/wvtransform/transformtospatialdomainwithf.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes
    + [`transformToSpatialDomainWithFAllDerivatives`](/classes/wvtransform/transformtospatialdomainwithfallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.
    + [`transformToSpatialDomainWithG`](/classes/wvtransform/transformtospatialdomainwithg.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes
    + [`transformToSpatialDomainWithGAllDerivatives`](/classes/wvtransform/transformtospatialdomainwithgallderivatives.html) transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.
    + [`transformUVEtaToWaveVortex`](/classes/wvtransform/transformuvetatowavevortex.html) transform fluid variables $$(u,v,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.
    + [`transformWaveVortexToUVWEta`](/classes/wvtransform/transformwavevortextouvweta.html) transform wave-vortex coefficients $$(A_+,A_-,A_0)$$ to fluid variables $$(u,v,\eta)$$.
+ Properties
  + [`effectiveHorizontalGridResolution`](/classes/wvtransform/effectivehorizontalgridresolution.html) returns the effective grid resolution in meters
+ Nonlinear flux and energy transfers
  + [`energyFluxFromNonlinearFlux`](/classes/wvtransform/energyfluxfromnonlinearflux.html) converts nonlinear flux into energy flux
  + [`nonlinearFlux`](/classes/wvtransform/nonlinearflux.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxForFlowComponents`](/classes/wvtransform/nonlinearfluxforflowcomponents.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxOperation`](/classes/wvtransform/nonlinearfluxoperation.html) The operation responsible for computing the nonlinear flux
  + [`nonlinearFluxWithGradientMasks`](/classes/wvtransform/nonlinearfluxwithgradientmasks.html) returns the flux of each coefficient as determined by the nonlinear flux operation
  + [`nonlinearFluxWithMask`](/classes/wvtransform/nonlinearfluxwithmask.html) returns the flux of each coefficient as determined by the nonlinear flux
+ Index Gymnastics
  + [`indexFromModeNumber`](/classes/wvtransform/indexfrommodenumber.html) return the linear index into a spectral matrix given (k,l,j)
  + [`isValidConjugateModeNumber`](/classes/wvtransform/isvalidconjugatemodenumber.html) returns a boolean indicating whether (k,l,j) is a valid conjugate mode number
  + [`isValidModeNumber`](/classes/wvtransform/isvalidmodenumber.html) returns a boolean indicating whether (k,l,j) is a valid mode number
  + [`isValidPrimaryModeNumber`](/classes/wvtransform/isvalidprimarymodenumber.html) returns a boolean indicating whether (k,l,j) is a valid primary (non-conjugate) mode number
+ Other
  + [`A0`](/classes/wvtransform/a0.html) geostrophic coefficients at reference time t0 (m)
  + [`A0Z`](/classes/wvtransform/a0z.html) 
  + [`Am`](/classes/wvtransform/am.html) negative wave coefficients at reference time t0 (m/s)
  + [`Ap`](/classes/wvtransform/ap.html) positive wave coefficients at reference time t0 (m/s)
  + [`ApmD`](/classes/wvtransform/apmd.html) 
  + [`ApmN`](/classes/wvtransform/apmn.html) 
  + [`K2`](/classes/wvtransform/k2.html) 
  + [`Lr2`](/classes/wvtransform/lr2.html) squared Rossby radius
  + [`PA0`](/classes/wvtransform/pa0.html) 
  + [`conjugateDimension`](/classes/wvtransform/conjugatedimension.html) 
  + [`dftBuffer`](/classes/wvtransform/dftbuffer.html) 
  + [`dftConjugateIndex`](/classes/wvtransform/dftconjugateindex.html) 
  + [`dftPrimaryIndex`](/classes/wvtransform/dftprimaryindex.html) 
  + [`diffX`](/classes/wvtransform/diffx.html) 
  + [`diffY`](/classes/wvtransform/diffy.html) 
  + [`diffZF`](/classes/wvtransform/diffzf.html) differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)
  + [`diffZG`](/classes/wvtransform/diffzg.html) differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)
  + [`dk`](/classes/wvtransform/dk.html) 
  + [`dl`](/classes/wvtransform/dl.html) 
  + [`dynamicalVariable`](/classes/wvtransform/dynamicalvariable.html) 
  + [`enstrophyFluxFromF0`](/classes/wvtransform/enstrophyfluxfromf0.html) 
  + [`h_0`](/classes/wvtransform/h_0.html) [Nj Nkl]
  + [`h_pm`](/classes/wvtransform/h_pm.html) [Nj Nkl]
  + [`hasMeanPressureDifference`](/classes/wvtransform/hasmeanpressuredifference.html) checks if there is a non-zero mean pressure difference between the top and bottom of the fluid
  + [`horizontalModes`](/classes/wvtransform/horizontalmodes.html) 
  + [`iOmega`](/classes/wvtransform/iomega.html) 
  + [`initializePrimaryFlowComponents`](/classes/wvtransform/initializeprimaryflowcomponents.html) 
  + [`isHydrostatic`](/classes/wvtransform/ishydrostatic.html) 
  + [`isequal`](/classes/wvtransform/isequal.html) 
  + [`kAxis`](/classes/wvtransform/kaxis.html) k coordinate
  + [`kl`](/classes/wvtransform/kl.html) dimension of the interleaved k-l wavenumber coordinate
  + [`lAxis`](/classes/wvtransform/laxis.html) l coordinate
  + [`modeNumberFromIndex`](/classes/wvtransform/modenumberfromindex.html) 
  + [`operationForDynamicalVariable`](/classes/wvtransform/operationfordynamicalvariable.html) This function is designed with the following goals:
  + [`qgpvFluxFromF0`](/classes/wvtransform/qgpvfluxfromf0.html) 
  + [`spectralVanishingViscosityFilter`](/classes/wvtransform/spectralvanishingviscosityfilter.html) Builds the spectral vanishing viscosity operator
  + [`summarizeDegreesOfFreedom`](/classes/wvtransform/summarizedegreesoffreedom.html) 
  + [`summarizeDynamicalVariables`](/classes/wvtransform/summarizedynamicalvariables.html) 
  + [`summarizeFlowComponents`](/classes/wvtransform/summarizeflowcomponents.html) 
  + [`totalEnergy`](/classes/wvtransform/totalenergy.html) 
  + [`totalEnergyOfFlowComponent`](/classes/wvtransform/totalenergyofflowcomponent.html) 
  + [`totalEnergySpatiallyIntegrated`](/classes/wvtransform/totalenergyspatiallyintegrated.html) 
  + [`totalEnstrophy`](/classes/wvtransform/totalenstrophy.html) 
  + [`totalEnstrophySpatiallyIntegrated`](/classes/wvtransform/totalenstrophyspatiallyintegrated.html) 
  + [`totalHydrostaticEnergy`](/classes/wvtransform/totalhydrostaticenergy.html) 
  + [`transformFromSpatialDomainWithFourier`](/classes/wvtransform/transformfromspatialdomainwithfourier.html) 
  + [`transformToSpatialDomainWithFourier`](/classes/wvtransform/transformtospatialdomainwithfourier.html) 
  + [`transformWithG_wg`](/classes/wvtransform/transformwithg_wg.html) 
  + [`variables`](/classes/wvtransform/variables.html) access the dynamical variables
  + [`variablesAtPosition`](/classes/wvtransform/variablesatposition.html) access the dynamical variables at any position in the domain
  + [`velocityField`](/classes/wvtransform/velocityfield.html) Return the velocity field, which is the sum of the gridded
  + [`version`](/classes/wvtransform/version.html) 
  + [`waveCoefficientsAtTimeT`](/classes/wvtransform/wavecoefficientsattimet.html) 
  + [`waveVortexTransformWithResolution`](/classes/wvtransform/wavevortextransformwithresolution.html) 
  + [`wvBuffer`](/classes/wvtransform/wvbuffer.html) 
  + [`wvConjugateIndex`](/classes/wvtransform/wvconjugateindex.html) 


---