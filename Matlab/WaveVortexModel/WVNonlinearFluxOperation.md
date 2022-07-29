Computes the nonlinear flux for a WVTransform

A WVNonlinearFluxOperation defines how energy moves between the wave-vortex coefficients---these are the nonlinear terms in the equations of motion, transformed into wave-vortex space. The most basic implementation is a freely evolving, unforced (and undamped) flux. Subclasses of WVNonlinearFluxOperation can implement custom forcing and damping.

A WVTransform is always initialized with a default nonlinear flux operation, but it can be overridden with,

```matlab
wvt.nonlinearFluxOperation = SomeNewNonlinearFluxOperation();
```

It is very likely you will want to use a custom nonlinear flux operation when integrating a model. In that case you would call,

```matlab
model = WVModel(wvt,nonlinearFlux=SomeNewNonlinearFluxOperation())
```

When creating a subclass of WVNonlinearFluxOperation, there are several important notes:

+ The output variables *must* be at least one of {Fp,Fm,F0}, in that order. The properties `doesFluxAp` etc. should be appropriately set to match the output.
+ You may also optionally output additional variables that are computed as a by product of your flux calculation. Those variables will then be cached, and will not have to be recomputed when needed.

- Declaration: classdef WVNonlinearFluxOperation < [WVOperation](/classes/wvoperation.html)
