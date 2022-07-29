---
layout: default
title: WVNonlinearFluxOperation
parent: Classes
has_children: false
has_toc: false
---

#  WVNonlinearFluxOperation

Computes the nonlinear flux for a WVTransform


---

## Declaration
```matlab
 classdef WVNonlinearFluxOperation < [WVOperation](/classes/wvoperation.html)
```

## Overview
  
  A WVNonlinearFluxOperation defines how energy moves between the
  wave-vortex coefficients---these are the nonlinear terms in the equations
  of motion, transformed into wave-vortex space. The most basic
  implementation is a freely evolving, unforced (and undamped) flux.
  Subclasses of WVNonlinearFluxOperation can implement custom forcing and
  damping.
  
  A WVTransform is always initialized with a default nonlinear flux
  operation, but it can be overridden with,
  
  ```matlab
  wvt.nonlinearFluxOperation = SomeNewNonlinearFluxOperation();
  ```
  
  It is very likely you will want to use a custom nonlinear flux operation
  when integrating a model. In that case you would call,
  
  ```matlab model =
  WVModel(wvt,nonlinearFlux=SomeNewNonlinearFluxOperation())
  ```
  
  When creating a subclass of WVNonlinearFluxOperation, there are several
  important notes:
  
  + The output variables *must* be at least one of {Fp,Fm,F0}, in that
  order. The properties `doesFluxAp` etc. should be appropriately set to
  match the output. + You may also optionally output additional variables
  that are computed as a by product of your flux calculation. Those
  variables will then be cached, and will not have to be recomputed when
  needed.
  
  


## Topics
+ Computation
  + [`compute`](/classes/wvnonlinearfluxoperation/compute.html) compute the promised variable
+ Equality
  + [`isequal`](/classes/wvnonlinearfluxoperation/isequal.html) check for equality with another nonlinear flux operation
+ Initialization
  + [`nonlinearFluxWithDoubleResolution`](/classes/wvnonlinearfluxoperation/nonlinearfluxwithdoubleresolution.html) create a new nonlinear flux operation with double the resolution
  + [`WVNonlinearFluxOperation`](/classes/wvnonlinearfluxoperation/wvnonlinearfluxoperation.html) create a new nonlinear flux operation
  + [`nonlinearFluxFromFile`](/classes/wvnonlinearfluxoperation/nonlinearfluxfromfile.html) initialize a nonlinear flux operation from NetCDF file
+ Properties
  + [`doesFluxAp`](/classes/wvnonlinearfluxoperation/doesfluxap.html) boolean indicating whether or not this operation returns Fp
  + [`doesFluxAm`](/classes/wvnonlinearfluxoperation/doesfluxam.html) boolean indicating whether or not this operation returns Fm
  + [`doesFluxA0`](/classes/wvnonlinearfluxoperation/doesfluxa0.html) boolean indicating whether or not this operation returns F0
  + [`name`](/classes/wvnonlinearfluxoperation/name.html) name of the operation
  + [`outputVariables`](/classes/wvnonlinearfluxoperation/outputvariables.html) array of WVVariableAnnotations describing the outputs of the computation
  + [`nVarOut`](/classes/wvnonlinearfluxoperation/nvarout.html) number of variables returned by the computation
  + [`f`](/classes/wvnonlinearfluxoperation/f.html) function handle to be called when computing the operation
+ Write to file
  + [`writeToFile`](/classes/wvnonlinearfluxoperation/writetofile.html) write information about the nonlinear flux operation to file


---