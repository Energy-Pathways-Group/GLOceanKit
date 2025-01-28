---
layout: default
title: Using the WVTransform
parent: Users guide
mathjax: true
nav_order: 2
has_toc: true
---

#  Using the WVTransform

At the time of this writing there are three `WVTransform` subclasses with different capabilities,
- `WVTransformConstantStratification` Non-hydrostatic, constant stratification
- `WVTransformHydrostatic` Hydrostatic, variable stratification
- `WVTransformSingleMode` Single mode, equivalent barotropic

These three `WVTransform` subclasses are used in the same way, but have three different requirements for initialization.

## Initialization

Constant stratification flows require you specify the domain size, number of grid points, and, optionally, the buoyancy frequency $$N_0$$ and latitude,
```matlab
wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0=N0,latitude=latitude);
```

The hydrostatic transformation requires you pass a function handle describing the stratification, e.g.,
```matlab
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=N2,latitude=latitude);
```

Note that you do not specify the grids, only the dimensions, as the grids are determined by the transforms. The $$x$$ and $$y$$ grids will always be evenly spaced grids appropriate for Fourier transforms, while the $$z$$ grid will be evenly spaced for constant stratification (sine and cosine transforms), but will have variable spacing when using variable stratification.

The equivalent barotropic transform requires that you specify an equivent depth $$h$$,
```matlab
wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.80,latitude=latitude);
```
where a typical oceanic value for the first barolinic mode would be around 80 cm.

## Initial conditions

Once a model is initialized, it is often the case that one would like to add initial conditions. The [WVTransform](/classes/wvtransform.html) class has numerous methods for adding initial conditions---including very generael initialization from any fluid state, as well as initialization specific to waves, inertial oscillations, and geostrophic motions.

### Initializing from $$(u,v,\eta)$$

The most direct methods for initializing the model are [`initWithUVEta`](/classes/wvtransform/initwithuveta.html) and [`initWithUVRho`](/classes/wvtransform/initwithuvrho.html) which take either $$(u,v,\eta)$$ or $$(u,v,\rho)$$. As a simple example, let's intialize with an inertial oscillation initial condition, $$(u_0 \exp(z/100),0,0)$$. In code, this is
```matlab
[X,Y,Z] = wvt.xyzGrid;
wvt.initWithUVEta( 0.2*exp(Z/100), 0*X, 0*X );
```

The call to `[X,Y,Z] = wvt.xyzGrid` returns three matrices, typically of size `[Nx Ny Nz]`, that contain the grid values. Behind the scenes this is calling the Matlab function `ndgrid`, but please do not create the grids yourself---a `WVTransform` may choose to order the the grid points differently than you expect.

These methods can be used to initialize from *any* flow fields that use the same stratification and boundary conditions as the `WVTransform` that is being used. For example, you might use output from another model and use [`initWithUVRho`](/classes/wvtransform/initwithuvrho.html) to initialize the WaveVortexMode with that output.

### Initializing waves, inertial oscillations, and geostrophic motions

The wave-vortex model provides methods for initializing the fluid with specific dynamical solutions, including inertial oscillations, internal gravity waves, and geostrophic (vortex) flows. 

To initialize individual waves, use  [`initWithWaveModes`](/classes/wvtransform/initwithwavemodes.html) and the related methods, e.g.,
```matlab
U = .2; phi=0;
omega = wvt.initWithWaveModes(10,0,1,phi,U,1);
period = 2*pi/omega;
```
will initialize a first baroclinic mode wave with a wavenumber of $$k_x = 10(2\pi)/L_x$$.

To initialize with a geostrophic stream function, use  [`initWithGeostrophicStreamfunction`](/classes/wvtransform/initwithgeostrophicstreamfunction.html) and the related methods, e.g.,
```matlab
Le = 35e3;
z0 = -wvt.Lz/4;
He = wvt.Lz/10;
U = 0.25; % m/s
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2) .* (erf((z-z0)/He)+1)/2;

wvt.setGeostrophicStreamfunction(psi);
```
creates a deep eddy.

The initialization methods for the [WVTransform](/classes/wvtransform.html) all use the same nomenclature,

- `init`---clears ALL variables `Ap`, `Am`, `A0`, then sets/adds
- `set`---clears only the component requested, and sets with new value.
- `add`---adds to existing component
- `removeAll`---remove all features of given type  

### Other initialization methods

It can also be useful to [`initWithRandomFlow`](/classes/wvtransform/initwithrandomflow.html) and [`removeEnergyFromAliasedModes`](/classes/wvtransform/removeenergyfromaliasedmodes.html) before doing a model run. Also check out the section on [reading and writing to file](/users-guide/reading-and-writing-to-file.html) for how to use [`initFromNetCDFFile`](/classes/wvtransform/initfromnetcdffile.html) to quickly read in the ocean state from a saved file.

## Examining the fluid state

Once you have a `WVTransform` instance, you can now query it to return different state variables. The term "state" is used because these variables tell you something about the state of the fluid at the time in question.

Standard variables that are available are documented with the [WVTransform](/classes/wvtransform.html) and include $$(u,v,w,\rho,p)$$, the sea-surface height and velocities, the quasi-geostrophic potential vorticity (qgpv), the nonlinear fluxes $$(F_p, F_m, F_0)$$ and others. 

There are two ways to access state variables. In the first way, you can query the transform using Matlab's dot syntax, e.g.,
```matlab
u = wvt.u;
```
will return a matrix of the $$u$$ velocity component. This can then be plotted in the usual way, e.g.,
```matlab
figure
pcolor(wvt.y,wvt.z,squeeze(u(1,:,:)).'), shading interp
```
where the properties of the WVTransform `wvt.y` and `wvt.z` refer to those dimensions.

The other way to access variables is to query them in a list using the [`variables`](/classes/wvtransform/variables.html) method,
 ```matlab
[u,v,w] = wvt.variableWithName('u','v','w');
```
this returns all the variables requested at the grid points. However, you can also the [`variableAtPositionWithName`](/classes/wvtransform/variablesatposition.html) method to return the value of any variable at *any* set of points in the domain. This method is used to do particle advection by the `WVModel`.


