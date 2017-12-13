InternalWaveModel
==============

The InternalWaveModel class creates linear internal wave model for user specified stratification.

If you use these classes, please cite the following paper,
- J. Early. Linear internal wave model, Ocean Modelling, in prep

### Table of contents
1. [Quick Start](#quick-start)
2. [Dynamical variable accessor methods](#dynamical-variable-accessor-methods)
3. [Built-in density profiles](#built-in-density-profiles)

------------------------

Quick start
------------

To initialize the linear internal wave model, you first have to set the problem dimensions,
```matlab
aspectRatio = 8;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

Nx = 64;
Ny = aspectRatio*64;
Nz = 65; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
```
and then initialization either the constant stratification model with some buoyancy frequency,
```matlab
N0 = 5.2e-3;
wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
```
or initialize the arbitrary stratification model with some density profile,
```matlab
[rho, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
z = linspace(min(zIn),max(zIn),100)';
maxModes = 32;
wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, maxModes, latitude);
```
The constant stratification model *must* model the entire ocean depth, because it uses the fast cosine transform to sum the vertical modes. The arbitrary stratification model, however, can be given `z` that only includes part of the domain, because it uses the slow (matrix multiplication) transform to sum the vertical modes. The argument `maxModes` sets how many vertical modes are used. The stratification profile in the example came from the [`InternalModes` class](../InternalModes/).

Once the model has been created, you can initialize the model with either a single plane wave,
```matlab
k0 = 4; % k=0..Nx/2
l0 = 0; % l=0..Ny/2
j0 = 20; % j=1..nModes, where 1 indicates the 1st baroclinic mode
U = 0.01; % m/s
sign = 1;

period = wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);
```
or a full spectrum of waves,
```matlab
wavemodel.InitializeWithGMSpectrum(1.0);
```

To access the velocity field or density field, you can simply request their value at any time,
```matlab
t = 0.0;
[u,v,w] = wavemodel.VelocityFieldAtTime(t);
rho = wavemodel.DensityFieldAtTime(t);
zeta = wavemodel.IsopycnalDisplacementFieldAtTime(t);
```

The model also supports external/free wave modes that do not fit in the gridded domain, as well as a series of functions for advecting particles.

Initialization
------------




Dynamical variable accessor methods
------------
Once the `InternalWaveModel` object is initialized there are a few useful accessor methods for access the velocity and density at different times and positions. The API uses the following terminology:

- *Eulerian* methods return values on the build in grid. They will contain the word `Field` in the name, such as `VelocityFieldAtTime` and will always return variables in arrays that match the grid dimensions.
- *Lagrangian* methods return values at arbitrary, user-specified positions. They will contain the phrase `Position` in the name, such as `VelocityAtTimePosition`.
- *Internal*, or *gridded*, refers to the portion of the wave field defined on the gridded field. For example,  `InternalVelocityFieldAtTime`.
- *External* refers to external wave modes; for example,  `ExternalVelocityFieldAtTime`.

The internal and external variables always sum to the total, which is what is returned by default.

The density is decomposed into two parts such that density = mean + perturbation. The mean is strictly a function of z and does not vary in time.

### Eulerian

The following methods return values on the grid points.

```matlab
t = 0.0;
[u,v,w] = wavemodel.VelocityFieldAtTime(t);
rho = wavemodel.DensityFieldAtTime(t);
zeta = wavemodel.IsopycnalDisplacementFieldAtTime(t);
```

The density field is defined such that each part can be accessed separately,
```matlab
rho_bar = wavemodel.DensityMeanField;
rho_prime = wavemodel.DensityPerturbationFieldAtTime(t);
```

Fundamentally these routines are just calling `VariableFieldsAtTime`, which lets you request any of the dynamical variables you want at a given time. If you are going for speed, it is best to call this function once for a given time point. You could request all variables like this,

```matlab
[u,v,w,rho_prime,zeta] = wavemodel.VariablesFieldsAtTime(t,'u','v','w','rho_prime','zeta');
```

### Lagrangian

These methods will return values at any point in space and time. The argument `interpolationMethod` specifies how off-grid values should be interpolated. Use `'exact'` for the slow, but accurate, spectral interpolation. Otherwise use `'spline'` or some other method accepted by Matlab's `interp` function.

As an example, lets create some floats at various depths in the water column.
```matlab
dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
nLevels = 5;
N = floor(wavemodel.Nx/3);
x = (0:N-1)*dx;
y = (0:N-1)*dy;
z = (0:nLevels-1)*(-wavemodel.Lz/(2*(nLevels-1)));
```
Now we can call the various Lagrangian methods to return the dynamical fields at those positions. First, we need to specify how we want to interpolate the gridded dynamical fields. The two most obvious choices are `'exact'`, if we want to use spectral interpolation, or `'spline'`, if we want a good accuracy speed tradeoff. It may also being worth using `'linear'` interpolation if accuracy is less important.

```matlab
interpolationMethod = 'spline';
[u,v,w] = wavemodel.VelocityAtTimePosition(t,x,y,z,interpolationMethod);
rho = wavemodel.DensityAtTimePosition(t,x,y,z,interpolationMethod);
zeta = wavemodel.IsopycnalDisplacementAtTimePosition(t,x,y,z,interpolationMethod);
```
As with the Eulerian accessor methods, you can request all the dynamical variables at the same time, e.g.
```matlab
[u,v,w,rho_prime,zeta] = wavemodel.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'u','v','w','rho_prime','zeta');
```
which provides the optimal speed advantage.

### Internal and external variables

If you so choose, you can also access the dynamical variables associated with the internal and external wave modes separately. The following code accesses the Eulerian and Lagrangian versions of the internal wave modes,
```matlab
[u,v,w,rho_prime,zeta] = wavemodel.InternalVariableFieldsAtTime(t,'u','v','w','rho_prime','zeta');
[u,v,w,rho_prime,zeta] = wavemodel.InternalVariablesAtTimePosition(t,x,y,z,interpolationMethod,'u','v','w','rho_prime','zeta');
```
In practice, the Eulerian function `InternalVariableFieldsAtTime` is doing all the heavy lifting by computing the FFTs for the gridded modes. The Lagrangian function `InternalVariablesAtTimePosition` uses the Eulerian solutions and interpolates the position at the user requested point, unless `'exact'` interpolation is requested.

The external dynamical variables can be access with,
```matlab
[u,v,w,rho_prime,zeta] = wavemodel.ExternalVariableFieldsAtTime(t,'u','v','w','rho_prime','zeta');
[u,v,w,rho_prime,zeta] = wavemodel.ExternalVariablesAtTimePosition(t,x,y,z,'u','v','w','rho_prime','zeta');
```
Notice that interpolation is not an option, because these values are always exact.

Built-in density profiles
------------

For testing purposes and convenience, there are a number of pre-defined density profiles that you can access. Simply call,
```matlab
[rho,N2,zIn] = InternalModes.StratificationProfileWithName(stratification)
```
where the variable `stratification` is a string. The returned values `rho` and `zIn` can be used directly as the first two arguments in initialization. The following options are available:
 - `'constant'` A constant stratification profile.
 - `'exponential'` The standard Garrett-Munk exponential profile.
 - `'pycnocline-constant'` Constant stratification with a pycnocline, taken from Cushman-Roisin & Beckers.
 - `'pycnocline-exponential'` Exponential profile with a deep pycnocline.
 - `'latmix-site1'` Attempts to recreate the full stratification profile (down to 5000 meters) at Latmix site 1, with an intense surface mixed layer and deep pycnocline.
 - `'latmix-site1-surface'` Recreates just the near surface features at Latmix Site 1.
 - `'latmix-site1-exponential'` Recreates the surface mixed layer and associated pycnocline, and then decays exponentially below 190m.
  - `'latmix-site1-constant'` Recreates the surface mixed layer and associated pycnocline, but then goes constant below 300m.
  

