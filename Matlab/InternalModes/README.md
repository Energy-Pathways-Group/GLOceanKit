InternalModes
==============

The InternalModes class can be used to quickly and accurately compute the vertical modes from arbitrary stratification.

The complete class hierarchy contains several implementations that include spectral methods, finite differencing, WKB approximated solutions, as well as the analytical solutions for constant stratification. The details are documented in Early, Lelong, and Smith (2018).

The classes contain many options, but also try to remain simple to use. Try `help InternalModes` in Matlab for a complete description, or use the Quick Start below.

Quick Start
------------

There are two primary methods of initializing this class: either
you specify density with gridded data, or you specify density as an
analytical function.

```matlab
im = InternalModes(rho,z,zOut,latitude);
```
'rho' must be a vector of gridded density data with length matching
'z'. zOut is the output grid upon which all returned function will
be given.

```matlab
modes = InternalModes(rho,zDomain,zOut,latitude);
```
'rho' must be a function handle, e.g.
```matlab
   rho = @(z) -(N0*N0*rho0/g)*z + rho0
   ```
zDomain must be an array with two values: `z_domain = [z_bottom z_surface];`

Once initialized, you can request variations of the density, e.g.,
```matlab
   N2 = modes.N2;
   rho_zz = modes.rho_zz;
   ```
or you can request the internal modes at a given wavenumber, k,
where k is 2*pi/wavelength.
```matlab
   [F,G,h] = modes.ModesAtWavenumber(0.01);
   ```
or frequency,
```matlab
   [F,G,h] = modes.ModesAtWavenumber(5*modes.f0);
   ```