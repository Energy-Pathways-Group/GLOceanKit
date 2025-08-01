Perform an operation and return a variable using a WVTransform

A WVOperation allows you to add new functionality to the WVTransform by defining one or more variables and creating an operation for computing those variables.

If the operation is simple and can be computed in one line, you can directly instantiate the WVOperation class and pass a function handle. For example,

```matlab
outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
wvt.addOperation(WVOperation('zeta_z',outputVar,f));
```

will enable direct calls to `wvt.zeta_z` to compute the vertical vorticity.

Note that a `WVOperation` that computes a single variable, must have the same name as the variable, as specified in `WVVariableAnnotation`. 

More involved calculations may require subclassing WVOperation and overriding the `compute` method. Note that every variable returned by the compute operation must be described with a `WVVariableAnnotation`.

- Declaration: classdef WVOperation < handle
