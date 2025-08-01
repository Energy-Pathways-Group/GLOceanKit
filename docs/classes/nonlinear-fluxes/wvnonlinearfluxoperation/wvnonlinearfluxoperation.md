---
layout: default
title: WVNonlinearFluxOperation
parent: WVNonlinearFluxOperation
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  WVNonlinearFluxOperation

create a new nonlinear flux operation


---

## Declaration
```matlab
 nlFluxOp = WVNonlinearFluxOperation(name,outputVariables,options)
```
## Parameters
+ `name`  name of the nonlinear flux operation
+ `outputVariables`  ordered list WVVariableAnnotation instances describing each variable returned by the compute method
+ `f`  (optional) directly pass a function handle, rather than override the compute method

## Returns
+ `nlFluxOp`  a new instance of WVNonlinearFluxOperation

## Discussion

  This class is intended to be subclassed, so it generally
  assumed that this initialization will not be called directly.
 
            
