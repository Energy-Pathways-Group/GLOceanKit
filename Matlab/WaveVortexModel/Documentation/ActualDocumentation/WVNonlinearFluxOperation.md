

The output variables *must* be at least one of {Fp,Fm,F0}, in that order. The properties `doesFluxAp` etc. should be appropriately set to match the output.

You may also optionally output additional variables that are computed as a by product of your flux calculation. Those variables will then be cached, making other function calls much quicker.

- Declaration: classdef WVNonlinearFluxOperation < WVOperation
