InternalModes
==============

The InternalModes class can be used to quickly and accurately compute the vertical modes from arbitrary stratification.

The complete class hierarchy contains several implementations that include spectral methods, finite differencing, WKB approximated solutions, as well as the analytical solutions for constant stratification. The details are documented in Early, Lelong, and Smith (2018).

The classes contain many options, but also try to remain simple to use. Try `help InternalModes` in Matlab for a complete description, or use the Quick Start below.

Quick Start
------------

The preferred method for initializing the InternalModes class is to define density as a function, e.g.,
```matlab
N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
g = 9.81;
rho0 = 1025; % density at the surface
b = 1300;
rho = @(z) rho0*(1 + b*N0*N0/(2*g)*(1 - exp(2*z/b)));
```
and then use that function to initialize the class,
```matlab
L = -5000;
zOut = linspace(L,0,200)';
latitude = 33;
im = InternalModes(rho,[L 0],zOut,latitude);
```
The `InternalModes` class was initialized with four arguments: the density function, an array specifying the domain (lower and upper boundary), the the output grid upon which all returned function will be given, and the latitude.

Now that the `im` object is initialized, you can request the internal modes at a given wavenumber, k, where k is 2*pi/wavelength.
```matlab
   [F,G,h,omega] = im.ModesAtWavenumber(2*pi/1000);
   ```
or frequency,
```matlab
   [F,G,h,k] = im.ModesAtWavenumber(5*modes.f0);
   ```
The arrays `F` and `G` contain the vertical modes for U/V and W, respectively. The arrays have dimensions `size(F)=[length(zOut) length(h)]`, meaning that each column `i` is a normal mode, `F(:,i)` with corresponding eigendepth `h(i)`.

Note that you can also request variations of the density, e.g.,
```matlab
   N2 = im.N2;
   rho_zz = im.rho_zz;
   ```
which will be on the same output grid `zOut` that you specified.


Although it is best to use initialize using a density function, you can also initialize using gridded data if needed. The syntax is nearly the same,
```matlab
N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
g = 9.81;
rho0 = 1025; % density at the surface
b = 1300;
z = linspace(-5000,0,500)';
rho = rho0*(1 + b*N0*N0/(2*g)*(1 - exp(2*z/b)));

im = InternalModes(rho,z,zOut,latitude);
```
where now you pass the gridded density data, `rho`, and its coordinate, `z`.

Convenience functions
------------
Once the `InternalModes` objects is initialized there are a few notable convenience functions. You can use,
```matlab
im.ShowLowestModesAtWavenumber(2*pi/1000)
```
and

```matlab
im.ShowLowestModesAtFrequency(5*modes.f0)
```
to quickly visualize the four lowest modes.


Built-in Density Profiles
------------

For testing purposes and convenience, there are a number of pre-defined density profiles that you can access. Simply call,
```matlab
[rho,N2,zIn] = InternalModes.StratificationProfileWithName(stratification)
```
where the variable `stratification` is a string. The returned values `rho` and `zIn` can be used directly as the first two arguments in initialization. The following options are available:
 - `'constant'` A constant stratification profile.
 - `'exponential'` The standard Garrett-Munk exponential profile.
 - `'pycnocline-constant'` Constant stratification with a pycnocline.
 - `'pycnocline-exponential'` Exponential profile with a deep pycnocline.
 - `'latmix-site1'` Attempts to recreate the full stratification profile (down to 5000 meters) at Latmix site 1, with an intense surface mixed layer and deep pycnocline.
 - `'latmix-site1-surface'` Just recreates the near surface features at Latmix Site 1.
 - `'latmix-site1-exponential'` Recreates the surface mixed layer and associated pycnocline, and then decays exponentially below 190m.
  - `'latmix-site1-constant'` Recreates the surface mixed layer and associated pycnocline, but then goes constant below 300m.
  

Advanced Features
----------------

Before delving into the more advanced features, you need to understand the class hierarchy used by InternalModes.

