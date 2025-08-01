---
layout: default
title: Creating new state variables
parent: Users guide
mathjax: true
nav_order: 4
---

#  Creating new state variables

The WaveVortexModel allows you to quickly and easily create custom operations to define new variables describing the state of the ocean. These new variables can saved to file and their values can be tracked along particle trajectories using the model.

## WVOperation

The simplest way to add new functionality to the WaveVortexModel is to create a new `WVOperation`. A `WVOperation` instance has two jobs: 1) describe the variable (or variables) that it is capable of producing and 2) implement an operation for computing those variables.

### Creating a new operation

As an example, lets add a new operation to compute relative vorticity. Assuming that `wvt` is a WVTransform instance,

```matlab
outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
wvt.addOperation(WVOperation('zeta_z',outputVar,f));
```

That's it! The `WVVariableAnnotation` class is used to describe the structure of the the output variable, and the function handle `f` performs the actual computation. The computation is not actually performed until it is needed.

Now that the WVTransform has a recipe for computing `zeta_z`, you can simply call `wvt.zeta_z` at any time. In fact, these operations compose---so you can create a new operation itself calls `wvt.zeta_z` as part of its computation.

### Saving to file

The new variable can be immediately saved to file because it already has a name, dimensions, and units.

Using the variable `zeta_z` from the example above, the transform can be written to file

```matlab
wvt.writeToFile('path/to/file.nc','zeta_z');
```

When doing a model simulation, this variable is also available to write to file as part of the time series see, e.g., `setNetCDFOutputVariables`.

### Tracking along particle trajectories

If the output variable has spatial dimension $$(x,y)$$ or $$(x,y,z)$$, then the WVModel can track its value along advected particle trajectories. See the documentation for `addParticles`, `setFloatPositions`, or `setDrifterPositions`. 

## Notes

In the example above a function handle was used to compute the relative vorticity because it only required a single line of code. However, it is often the case that computations will be more involved. In that case, you should subclass the `WVOperation`.

If a `WVOperation` produces a variable with the same name as one that already exists, it will *replace* the existing `WVOperation`. This is both useful and dangerous.

A `WVOperation` can return multiple variables.

All computed variables are cached, so as not to repeat unnecessary calculations.
