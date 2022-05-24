---
layout: default
title: WaveVortexTransform
parent: Classes
---
#  WaveVortexTransform

The WaveVortexTransform subclasses encapsulate data representing the state of the ocean at a given instant in time (e.g., u, v, w, and rho). What makes the WaveVortexTransform subclasses special is that the state of the ocean is represented as energetically independent waves and geostrophic motions (vortices). These classes can be queried for other information, e.g., Ertel PV, relative vorticity, etc.

## Overview



## Topics
+ [Initialization](#initialization)
  + [`WaveVortexTransformConstantStratification`](#dimensionwithname)
  + [`WaveVortexTransformHydrostatic`](#dimensionwithname)
  + [`InitFromFile`](#dimensionwithname)
  + [`waveVortexTransformWithResolution`](#dimensionwithname)
+ [Ocean State](#ocean-state)
  + [`Variables`](#attributes)
  + [`VariablesAtPosition`](#attributes)
  + [`VelocityField`](#attributes)
  + Ocean State Variables
    + [Ap, Am, A0](#attributes)
    + [u, v, w, eta, rho_prime](#attributes)
    + [`rho_total`](#attributes)
    + [Fp, Fm, F0](#attributes)
    + [`Fm`](#attributes)
    + [`F0`](#attributes)
  + [`addTransformOperation`](#addtransformoperation)
+ [Initial Conditions](#initial-conditions)
  + Waves
    + [`InitializeWithPlaneWave`](#attributes)
    + [`RemoveAllGriddedWaves`](#attributes)
    + [`SetGriddedWavesWithWavemodes`](#attributes)
    + [`AddGriddedWavesWithWavemodes`](#attributes)
    + [`WaveCoefficientsFromGriddedWaves`](#attributes)
    + [`InitializeWithGMSpectrum`](#attributes)
    + [`InitializeWithSpectralFunction`](#attributes)
  + Inertial Oscillations
    + [`SetInertialMotions`](#attributes)
  + Geostrophic Motions
    + [`SetGeostrophicStreamfunction`](#attributes)
+ [Energetics](#energetics)
  + [`totalEnergy`](#totalenergy)
  + [`totalSpectralEnergy`](#addattributenamevalue)
  + [`totalSpectralEnergy`](#addattributenamevalue)
  + Major Constituents
    + [`inertialEnergy`](#addattributenamevalue)
    + [`waveEnergy`](#addattributenamevalue)
    + [`geostrophicEnergy`](#addattributenamevalue)
  + Geostrophic Constituents
    + [`barotropicGeostrophicEnergy`](#addattributenamevalue)
    + [`baroclinicGeostrophicEnergy`](#addattributenamevalue)
  + Inertia-Gravity Wave Constituents
    + [`barotropicInertialEnergy`](#addattributenamevalue)
    + [`baroclinicInertialEnergy`](#addattributenamevalue)
    + [`internalWaveEnergyPlus`](#addattributenamevalue)
    + [`internalWaveEnergyMinus`](#addattributenamevalue)
  + [`summarizeEnergyContent`](#addattributenamevalue)
+ [Constituent Flow Masks](#energetics)
  + [`MasksForFlowConstituents`](#totalenergy)
  + [`MasksForAllFlowConstituents`](#addattributenamevalue)
  + [`MaskForAliasedModes`](#addattributenamevalue)
+ [Read And Write](#read-and-write)
  + [`WriteToFile`](#dimensions)
  + [`InitFromFile`](#dimensionwithname)

---

## Initialization

---

### `WaveVortexTransformConstantStratification(dims, n, latitude, N0, rho0, varargin)`
> Create a new...
>
> Usage


---

## Ocean State

---

### `Variables`
> A `containers.Map` type that contains the key-value pairs of all global attributes in the NetCDF file. This is intended to be *read only*. If you need to add a new attribute to file, use [`addAttribute`](#addattribute).
>
> Usage
> ```matlab
> model = ncfile.attributes('model');
> ```

---

### `VariablesAtPosition(name,value)`
> Adds a global attribute to the NetCDF file.
>
> Usage
> `ncfile.addAttribute('model','WaveVortexModel');`
>
> Parameters
> - `name` (key) string with the name of the property
> - `value` value




