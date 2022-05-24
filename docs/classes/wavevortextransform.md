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
  + [`InitFromFile`](#dimensionwithname)
  + [`waveVortexTransformWithResolution`](#dimensionwithname)
+ [Ocean State](#ocean-state)
  + [`Variables`](#attributes)
  + [`VariablesAtPosition`](#attributes)
  + [`VelocityField`](#attributes)
  + Ocean State Variables
    + [`Ap`](#attributes)
    + [`Am`](#attributes)
    + [`A0`](#attributes)
    + [`u`](#attributes)
    + [`v`](#attributes)
    + [`w`](#attributes)
    + [`eta`](#attributes)
    + [`rho_prime`](#attributes)
    + [`rho_total`](#attributes)
    + [`Fp`](#attributes)
    + [`Fm`](#attributes)
    + [`F0`](#attributes)
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
+ [Read And Write](#read-and-write)
  + [`WriteToFile`](#dimensions)
  + [`InitFromFile`](#dimensionwithname)

---

## Initialization

---

### `self = NetCDFFile(path,overwriteExisting)`
> Open an existing file or create a new empty file.
>
> Usage
> ```matlab
> ncfile = NetCDFFile(path)
> ```
> will load an existing file (if one exists) or create a new file (if none exists).
> ```matlab
> ncfile = NetCDFFile(path,'OVERWRITE_EXISTING')
> ```
> will delete any existing file and create a new file.

---

## Global Attributes

---

### `attributes`
> A `containers.Map` type that contains the key-value pairs of all global attributes in the NetCDF file. This is intended to be *read only*. If you need to add a new attribute to file, use [`addAttribute`](#addattribute).
>
> Usage
> ```matlab
> model = ncfile.attributes('model');
> ```

---

### `addAttribute(name,value)`
> Adds a global attribute to the NetCDF file.
>
> Usage
> `ncfile.addAttribute('model','WaveVortexModel');`
>
> Parameters
> - `name` (key) string with the name of the property
> - `value` value

---

## Coordinate Dimensions

---

### `dimensions`
> An array of NetCDFDimension objects for each coordinate dimension defined in the NetCDF file. The dimensions order in the array should reflect the underlying dimensionID defined in the NetCDF file.
>
> Usage
> ```matlab
> dim = ncfile.dimensions(dimID+1); % get the dimension with dimID
> ```

---

### `dimensionWithName`
> An container.Map of NetCDFDimension objects keyed by dimension name.
>
> Usage
> ```matlab
> xDim = ncfile.dimensionWithName('x');
> ```

---

### `[dimension,variable] = addDimension(name,value,properties,dimLength)`
> Adds a both a new dimension and its associated coordinate variable to the NetCDF file.
> 
> Usage
> ```matlab
> x = linspace(0,10,11);
> ncfile.addDimension('x',x,[]);
> ```
>
> Parameters
> - `name` string with the name of the dimension
> - `value` array of values along that dimension, or empty
> - `properties` containers.Map containing any key-value pairs to be associated with the dimension.
> - `dimLength` 
>
> Return Value
> - `dimension` a `NetCDFDimension` object with the newly create dimension
> - `variable` a `NetCDFVariable` object with the associated coordinate variable


---

### `[dimension,variable] = addMutableDimension(name,properties)`
> Add a new mutable dimension with coordinate variable to the NetCDF file.

---




