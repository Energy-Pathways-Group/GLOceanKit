---
layout: default
title: WaveVortexModel
parent: Classes
has_children: true
has_toc: false
---

#  WaveVortexModel

The WaveVortexModel is responsible for time-stepping (integrating) the ocean state forward in time using a WaveVortexTransform.

## Discussion
 
  model = WaveVortexModel(wvt) creates a new model using the wvt
  (WaveVortexTransform)
 
            


## Topics
+ [Initialization](#initialization)
  + [`modelFromFile`](/classes/wavevortexmodel/modelfromfile.html) Initialize a model from an existing file
  + [`WaveVortexModel`](/classes/wavevortexmodel/wavevortexmodel.html) Initialize a model from a WaveVortexTransform instance
+ [Model Properties](#model-properties)
  + [`wvt`](/classes/wavevortexmodel/wvt.html) The WaveVortexTransform instance the represents the ocean state.
  + [`nonlinearFlux`](/classes/wavevortexmodel/nonlinearflux.html) The operation responsible for computing the nonlinear flux of the model
  + [`linearDynamics`](/classes/wavevortexmodel/lineardynamics.html) Indicates whether or not the model is using linear or nonlinear dynamics.
  + [`t`](/classes/wavevortexmodel/t.html) Current model time (seconds)
+ [Integration](#integration)
  + [`integrateToNextOutputTime`](/classes/wavevortexmodel/integratetonextoutputtime.html) Time step the model forward to the next output time
  + [`integrateToTime`](/classes/wavevortexmodel/integratetotime.html) Time step the model forward to the requested time.
  + [`setupIntegrator`](/classes/wavevortexmodel/setupintegrator.html) Customize the time-stepping
  + [`outputInterval`](/classes/wavevortexmodel/outputinterval.html) Model output interval (seconds)
+ [Particles](#particles)
  + [`drifterPositions`](/classes/wavevortexmodel/drifterpositions.html) Current positions of the drifter particles
  + [`setDrifterPositions`](/classes/wavevortexmodel/setdrifterpositions.html) Set positions of drifter-like particles to be advected.
  + [`floatPositions`](/classes/wavevortexmodel/floatpositions.html) Returns the positions of the floats at the current time as well as the value of the fields being tracked.
  + [`setFloatPositions`](/classes/wavevortexmodel/setfloatpositions.html) Set positions of float-like particles to be advected by the model.
  + [`particlePositions`](/classes/wavevortexmodel/particlepositions.html) Positions and values of tracked fields of particles at the current model time.
  + [`addParticles`](/classes/wavevortexmodel/addparticles.html) Add particles to be advected by the flow.
+ [Tracer](#tracer)
  + [`tracer`](/classes/wavevortexmodel/tracer.html) Scalar field of the requested tracer at the current model time.
  + [`addTracer`](/classes/wavevortexmodel/addtracer.html) Add a scalar field tracer to be advected by the flow
+ [Writing to NetCDF files](#writing-to-netcdf-files)
  + [`createNetCDFFileForModelOutput`](/classes/wavevortexmodel/createnetcdffileformodeloutput.html) Create a NetCDF file for model output
  + [`removeNetCDFOutputVariables`](/classes/wavevortexmodel/removenetcdfoutputvariables.html) Remove variables from the list of variables to be written to the NetCDF variable during the model run.
  + [`setNetCDFOutputVariables`](/classes/wavevortexmodel/setnetcdfoutputvariables.html) Set list of variables to be written to the NetCDF variable during the model run.
  + [`addNetCDFOutputVariables`](/classes/wavevortexmodel/addnetcdfoutputvariables.html) Add variables to list of variables to be written to the NetCDF variable during the model run.
  + [`ncfile`](/classes/wavevortexmodel/ncfile.html) Reference to the NetCDFFile being used for model output
  + [`netCDFOutputVariables`](/classes/wavevortexmodel/netcdfoutputvariables.html) List of all StateVariables being written to NetCDF file


---