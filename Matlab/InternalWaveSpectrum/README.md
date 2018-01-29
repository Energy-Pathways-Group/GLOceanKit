InternalWaveSpectrum
==============

These class show the Garrett-Munk spectrum and its various approximations for different stratifications.

If you use these classes, please cite the following paper,
- J. Early, et al., Variance relations of the Garrett-Munk spectrum. Journal of Marine Research, in prep.

### Table of contents
1. [Initialization](#initialization)
2. [Variances](#variances)
3. [Spectra](#spectra)

------------------------


Initialization
------------

To initialize the `GarrettMunkSpectrum` class you must first pass a density profile in exactly the same forms as the [`InternalModes` class](../InternalModes/). For example, to initialize with an exponential profile, you would define the density function as follows,
```matlab
N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
g = 9.81; % gravity, in m/s^2
rho0 = 1025; % density at the surface
b = 1300; % decay scale, inmeters
L = 5000; % depth, in meters

rho = @(z) rho0*(1 + b*N0*N0/(2*g)*(1 - exp(2*z/b)));
```
and then use that function to initialize the class,
```matlab
latitude = 33;
GM = GarrettMunkSpectrum(rho,[-L 0],latitude);
```
Note that this class will compute the internal modes across the entire frequency spectrum, and so it may take several minutes to initialize. Once computed, these modes are stored in the `GM` object, so *do not overwrite* the `GM` object once these modes have been computed. Reuse the object.

Variances
------------
Now that the class is intialized, you can retrieve variances as a function of depth. For example,
```matlab
z = linspace(-L,0,5000)';

Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.IsopycnalVariance(z);
Ew = GM.VerticalVelocityVariance(z);
```
All values are in SI units.

Spectra
------------
You can also retrieve the distributions of the variances, i.e., the spectra of the these same quantities at a given depth. So let's grab the frequency spectra for these quantities at three different depths.
```matlab
depths = [0 -b/2 -b];
omega = linspace(-N0,N0,200);

Suv = GM.HorizontalVelocitySpectrumAtFrequencies(depths,omega);
Siso = GM.IsopycnalSpectrumAtFrequencies(depths,omega);
Sw = GM.VerticalVelocitySpectrumAtFrequencies(depths,omega);
```

Approximations
------------

By default, the variances and spectra are computed using the *exact* variance relations for the Garrett-Munk spectrum. However, it is often useful to consider the various approximations that are used. The three variance and three spectral functions above can all be given an optional argument specifying the approximation to be used, wit the following valid values,
- `'exact'` Uses the exact modal summation.
- `'wkb'` Uses wkb approximated vertical modes.
- `'wkb-hydrostatic'` Uses wkb approximated modes with the hydrostatic approximation, following Levine (2002).
- `'gm'` Uses the standard Garrett-Munk scalings, often referred to as the 'wkb scaled' variance relations.
