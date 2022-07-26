WaveVortexModel
==============

The WaveVortexModel decomposes a Boussinesq fluid with arbitrary stratification into geostrophic motions and interia-gravity waves. This model is complete, and therefore also can be time-stepped forward

This model has several significant features.

1. The internal waves can be initialized with Garrett-Munk spectrum, consistent for bounded domains and arbitrary stratification.
2. The model can be initialized from variables (u,v,rho), and thus performs a wave-vortex decomposition following the methodology in [Early, Lelong and Sundermeyer, 2021 JFM](https://arxiv.org/abs/2002.06267).

If you use these classes, please cite the JFM paper.

### Table of contents
1. [Quick Start](#quick-start)
2. [Introduction](#introduction)
4. [Gridded wave modes](#gridded-wave-modes)
5. [External wave modes](#external-wave-modes)
6. [Dynamical variables](#dynamical-variables)
    1. [Eulerian](#eulerian)
    2. [Lagrangian](#lagrangian)
    3. [Internal and external variables](#internal-and-external-variables)
7. [Initialization](#initialization)

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

period = wavemodel.setWaveModes(k0,l0,j0,U,sign);
```
or a full spectrum of waves,
```matlab
wavemodel.initWithGMSpectrum(1.0);
```

To access the velocity field or density field, you can simply request their value at any time,
```matlab
t = 0.0;
[u,v,w] = wavemodel.VelocityFieldAtTime(t);
rho = wavemodel.DensityFieldAtTime(t);
zeta = wavemodel.IsopycnalDisplacementFieldAtTime(t);
```

The model also supports external/free wave modes that do not fit in the gridded domain, as well as a series of functions for advecting particles.

Introduction
------------

The linear internal wave model adds together an arbitrary number single-mode plane-wave solutions to the linearized equations of motion on the f-plane. The model is entirely analytical---there is no numerical time-stepping routine---and the solutions are simply evaluated at the requested time.

The `InternalWaveModel` class is divided into two subclasses, depending on the stratification. The constant stratification model `InternalWaveModelConstantStratification` takes advantage of the fact that the vertical modes in constant stratification are simply sines and cosines, and uses the fast fourier transform to sum the vertical modes. For more general stratification profiles, the class `InternalWaveModelArbitraryStratification` computes the vertical modes using the [`InternalModes` class](../InternalModes/) and then sums them with a (slow) matrix multiplication. This technique will often be *slower* than simplying time-stepping the linearized equations of motion using a Runge-Kutta time stepping routine, but comes with the advantage that the vertical extent of the model can be arbitrarily defined. This allows a model to be constructed with high resolution in regions of interested.

Initializing the model creates a user-specified standard evenly-spaced, periodic grid. This grid useful because it allows wave modes to be computed with the Fast Fourier Transform (FFT) algorithm and both the constant and arbitrary stratification classes use this grid to quickly compute the horizontal component of the solution. However, such a grid limits the wavelengths of each wave to be integer multiples of the grid spacing. To get around this limitation, the linear wave model supports an arbitrary number of [external (or free) wave modes](#external-wave-modes) that are not required to fit in the gridded model. These wave modes are linearly added to the gridded modes, but must be computed using the slower discrete fourier transform.

Gridded wave modes
------------

The gridded wave modes by individually added or set by specifying the mode number and amplitude,
```matlab
[omega,k,l] = wavemodel.addWaveModes(kMode, lMode, jMode, phi, Amp, signs);
[omega,k,l] = wavemodel.setWaveModes(kMode, lMode, jMode, phi, Amp, signs);
```
where  `-Nx/2 < kMode < Nx/2` and `-Ny/2 < lMode < Ny/2` are integers specifying which modes you want to initialize. The `Set` method will remove all gridded modes and then add the list you give it, and the `Add` method will append these modes to the existing list.

Equivalently, you can call
```matlab
period = wavemodel.setWaveModes(k0,l0,j0,U,sign);
```
which just calls `setWaveModes` with a phase of 0.

The clear the gridded wave modes, call
```matlab
wavemodel.removeAllWaves();
```
and they will be set to 0 amplitude.


External wave modes
------------

The external wave modes can be added or set either by specifying the wavelength of the modes,
```matlab
omega = wavemodel.addExternalWavesWithWavenumbers(k,l,j,phi,A,norm);
omega = wavemodel.setExternalWavesWithWavenumbers(k,l,j,phi,A,norm);
```
or by specifying the frequency of the modes,
```matlab
k = wavemodel.addExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
k = wavemodel.setExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
```
The `Set` methods will remove all external modes and then add the list you give it, and the `Add` methods will append these modes to the existing list. You can call,
```matlab
wavemodel.removeAllExternalWaves();
```
to remove all external modes.

The amplitude of the waves is set with respect to a given norm of the internal mode (see [`InternalModes` class](../InternalModes/) for more details). Typically, one would want to use `Normalization.uMax`  if you are simply setting the maximum wave velocity.

You can extract the nonzero *gridded* wave components, and feed those directly into the external wave modes. The amplitude in this case uses `Normalization.kConstant`. This can be done with two lines of code,
```matlab
[omega, alpha, mode, phi, A, norm] = wavemodel.waveModesFromWaveCoefficients();
wavemodel.setExternalWavesWithFrequencies(omega, alpha, mode, phi, A, norm);
```
but of course now have you external waves with the exact same values as the gridded waves, which isn't likely very helpful.

Dynamical variables
------------
Once the `InternalWaveModel` object is initialized there are a few useful accessor methods for access the velocity and density at different times and positions. The API uses the following terminology:

- *Eulerian* methods return values on the built-in grid. They will contain the word `Field` in the name, such as `VelocityFieldAtTime` and will always return variables in arrays that match the grid dimensions.
- *Lagrangian* methods return values at arbitrary, user-specified positions. They will contain the phrase `Position` in the name, such as `VelocityAtTimePosition`.

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
Now we can call the various Lagrangian methods to return the dynamical fields at those positions. First, we need to specify how we want to interpolate the gridded dynamical fields. The two most obvious choices are `'exact'`, if we want to use spectral interpolation, or `'spline'`, if we want a good accuracy speed tradeoff. It may also being worth using `'linear'` interpolation if accuracy is less important. The following code returns the values of the dynamical variables at the positions we requested,

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

As explained in the introduction, the model is composed of *gridded* or *internal* wave modes and *external* or *free* wave modes.

If you so choose, you can also access the dynamical variables associated with the internal and external wave modes separately (which always sum to the total). The following code accesses the Eulerian and Lagrangian versions of the internal wave modes,
```matlab
[u,v,w,rho_prime,zeta] = wavemodel.InternalVariableFieldsAtTime(t,'u','v','w','rho_prime','zeta');
[u,v,w,rho_prime,zeta] = wavemodel.InternalVariablesAtTimePosition(t,x,y,z,interpolationMethod,'u','v','w','rho_prime','zeta');
```
In practice, the Eulerian function `InternalVariableFieldsAtTime` is doing all the heavy lifting by computing the FFTs for the gridded modes. The Lagrangian function `InternalVariablesAtTimePosition` uses the Eulerian solutions and interpolates the position at the user requested point, unless `'exact'` interpolation is requested.

The external dynamical variables can be access with,
```matlab
[u,v,w,rho_prime,zeta] = wavemodel.externalVariableFieldsAtTime(t,'u','v','w','rho_prime','zeta');
[u,v,w,rho_prime,zeta] = wavemodel.externalVariablesAtTimePosition(t,x,y,z,'u','v','w','rho_prime','zeta');
```
Notice that interpolation is not an option, because these values are always exact.

Initialization
------------

  The model can be initialized with a Garrett-Munk spectrum by calling,
  ```matlab
  wavemodel.initWithGMSpectrum(1.0);
  ```
where the only required argument indicates the GM reference level. This function also takes the following name/value pairs.

- `'j_star'` takes any integer value, default is 3.
- `'shouldRandomizeAmplitude'` takes a 0 or 1 to indicate whether or not the amplitude should be randomized with a Gaussian random variable with expectation matching the GM value. The default is 0.
- `'maxDeltaOmega'` is the maximum width in frequency that will be integrated over for assigned energy. By default it is self.Nmax-self.f.
- `'initializeModes'`  is used to determine which modes get initialized. Possible values are `'all'`, `'internalOnly'`, or `'externalOnly'` . Default is `'all'`
- `'energyWarningThreshold'` will provide a warning if the energy of a single mode exceeds a certain value of the total energy in that modal band. Values between 0 and 1. Default is 0.5 (e.g., you get a warning if the energy in a single mode exceeds 50% of the total energy).
- `'excludeNyquist'` takes a 0 or 1 to indicate whether or not to include the Nyquist wavenumbers in initialization. Default 1.
- `'maxK'` and `'minK'` can be set to limit the wavenumbers that will be initialized. By default all resolved wavenumbers are initialized (other than the Nyquist, as above).
- `'minMode'` and `'maxMode'` can be set to limit the vertical modes that will be initialized. By default all resolved vertical modes are initialized.
