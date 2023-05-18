---
layout: default
title: WVOperation
parent: WVOperation
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  WVOperation

create a new WVOperation for computing a new variable


---

## Declaration
```matlab
 operation = WVOperation(name,outputVariables,f)
```
## Parameters
+ `name`  name of the operation
+ `outputVariables`  ordered array of WVVariableAnnotations
+ `f`  function handle that takes a WVTransform as an argument and returns variables matching outputVariable.

## Discussion

            - Return operation: a new WV operation instance
