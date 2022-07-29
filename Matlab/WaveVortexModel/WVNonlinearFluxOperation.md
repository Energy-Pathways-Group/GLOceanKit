A WVNonlinearFluxOperation defines how energy moves between the wave-vortex coefficients---these are the nonlinear terms in the equations of motion, transformed into wave-vortex space. The most basic implementation is a freely evolving, unforced (and undamped) flux. Subclasses of WVNonlinearFluxOperation can implement custom forcing and damping.

The output variables *must* be at least one of {Fp,Fm,F0}, in that order. The properties `doesFluxAp` etc. should be appropriately set to match the output.

You may also optionally output additional variables that are computed as a by product of your flux calculation. Those variables will then be cached, making other function calls much quicker.

- Declaration: classdef WVNonlinearFluxOperation < WVOperation
