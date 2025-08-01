---
layout: default
title: Reading and writing to file
parent: Users guide
mathjax: true
nav_order: 3
has_toc: true
---

#  Reading and writing to file

The `WVTransform` and the `WVModel` both read and write to NetCDF files.

## WVTransform

After you create a `WVTransform` instance
```matlab
wvt = WVTransformConstantStratification([50e3 50e3 1300], [64 64 32]);
```
you call [`writeToFile`](/classes/wvtransform/writetofile.html)
```matlab
wvt.writeToFile('test.nc');
```
and all the information needed to exactly re-create the `wvt` instance will be written to file.

To create a new `WVTransform` instance from file, call the static method [`waveVortexTransformFromFile`](/classes/wvtransform/wavevortextransformfromfile.html)
```matlab
wvt2 = WVTransform.waveVortexTransformFromFile('test.nc');
```
The two instances `wvt` and `wvt2` are now equivalent.

### Adding variables

Any variable that the `WVTransform` instance knows about (including [custom variables](/users-guide/operations.html)) can also be written to file, but including its name in the variable argument list. For example, if you want to add the variables $$u$$, $$v$$, and $$\zeta_z$$ then call,
```matlab
wvt.writeToFile('test.nc','u','v','zeta_z');
```
and the data will be written, accessible from anything that reads NetCDF files.

### Reading model output

Model output contains multiple time points from which the a `WVTransform` instance can be initialized. To initialize a `WVTransform` instance from a specific point in the model output, call
```matlab
wvt = WVTransform.waveVortexTransformFromFile('test.nc',iTime=100);
``` 
where `iTime=100` indicates the index along the time dimension.

It is often the case that for analysis of model output, you want to read the model output at multiple time points. You *could* intialize a new `WVTransform` instance at each time index, but this involves unnecessary computation. Instead, you should update the existing `WVTransform` instance `wvt` with the data from that time point.

First load the NetCDF file, then initialize the existing instance from a specific time point in that file with [`initFromNetCDFFile`](/classes/wvtransform/initfromnetcdffile.html). For example,
```matlab
ncfile = NetCDFFile('test.nc');
wvt.initWithFile(ncfile,iTime=100);
```

## WVModel

Time series model output is created with

```matlab
model = WVModel(wvt,nonlinearFlux=WVNonlinearFlux(wvt,shouldAntialias=1));
model.setupIntegrator(timeStepConstraint="min",outputInterval=wvt.inertialPeriod/10);
model.createNetCDFFileForModelOutput("test.nc",shouldOverwriteExisting=0);
model.integrateToTime(10*wvt.inertialPeriod);
```

which can then be read in using the above commands.
