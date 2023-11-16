---
layout: default
title: WVModel
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 2
---

#  WVModel

The WVModel is responsible for time-stepping (integrating) the ocean state forward in time, as represented by a WVTransform.


---

## Overview
 
  Assuming you have already initialized a WVTransform, e.g.,
  ```matlab
  wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);
  ```
  and maybe set some initial conditions, you can then initialize the
  model,
  ```matlab
  model = WVModel(wvt)
  ```
  
  By default the model only takes a linear time-step. To specify a
  nonlinear flux on initialization, for example,
 ```matlab
  model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,u_damp=wvt.uMax));
 ```
 
  You can also initialize a model from existing output,
  ```matlab
  model = WVModel.modelFromFile('SomeFile.nc');
 ```
  
            


## Topics
+ Initialization
  + [`WVModel`](/classes/wvmodel/wvmodel.html) Initialize a model from a WVTransform instance
  + [`modelFromFile`](/classes/wvmodel/modelfromfile.html) Initialize a model from an existing file
+ Model Properties
  + [`linearDynamics`](/classes/wvmodel/lineardynamics.html) Indicates whether or not the model is using linear or nonlinear dynamics.
  + [`nonlinearFluxOperation`](/classes/wvmodel/nonlinearfluxoperation.html) The operation responsible for computing the nonlinear flux of the model
  + [`t`](/classes/wvmodel/t.html) Current model time (seconds)
  + [`wvt`](/classes/wvmodel/wvt.html) The WVTransform instance the represents the ocean state.
+ Integration
  + [`integrateToNextOutputTime`](/classes/wvmodel/integratetonextoutputtime.html) Time step the model forward to the next output time
  + [`integrateToTime`](/classes/wvmodel/integratetotime.html) Time step the model forward to the requested time.
  + [`outputIndex`](/classes/wvmodel/outputindex.html) output index of the current/most recent step.
  + [`outputInterval`](/classes/wvmodel/outputinterval.html) Model output interval (seconds)
  + [`setupIntegrator`](/classes/wvmodel/setupintegrator.html) Customize the time-stepping
+ Particles
  + [`addParticles`](/classes/wvmodel/addparticles.html) Add particles to be advected by the flow.
  + [`drifterPositions`](/classes/wvmodel/drifterpositions.html) Current positions of the drifter particles
  + [`floatPositions`](/classes/wvmodel/floatpositions.html) Returns the positions of the floats at the current time as well as the value of the fields being tracked.
  + [`particlePositions`](/classes/wvmodel/particlepositions.html) Positions and values of tracked fields of particles at the current model time.
  + [`setDrifterPositions`](/classes/wvmodel/setdrifterpositions.html) Set positions of drifter-like particles to be advected.
  + [`setFloatPositions`](/classes/wvmodel/setfloatpositions.html) Set positions of float-like particles to be advected by the model.
+ Tracer
  + [`addTracer`](/classes/wvmodel/addtracer.html) Add a scalar field tracer to be advected by the flow
  + [`tracer`](/classes/wvmodel/tracer.html) Scalar field of the requested tracer at the current model time.
+ Writing to NetCDF files
  + [`addNetCDFOutputVariables`](/classes/wvmodel/addnetcdfoutputvariables.html) Add variables to list of variables to be written to the NetCDF variable during the model run.
  + [`createNetCDFFileForModelOutput`](/classes/wvmodel/createnetcdffileformodeloutput.html) Create a NetCDF file for model output
  + [`ncfile`](/classes/wvmodel/ncfile.html) Reference to the NetCDFFile being used for model output
  + [`netCDFOutputVariables`](/classes/wvmodel/netcdfoutputvariables.html) List of all StateVariables being written to NetCDF file
  + [`removeNetCDFOutputVariables`](/classes/wvmodel/removenetcdfoutputvariables.html) Remove variables from the list of variables to be written to the NetCDF variable during the model run.
  + [`setNetCDFOutputVariables`](/classes/wvmodel/setnetcdfoutputvariables.html) Set list of variables to be written to the NetCDF variable during the model run.
+ Other
  + [`recordNetCDFFileHistory`](/classes/wvmodel/recordnetcdffilehistory.html) 


---