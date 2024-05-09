---
layout: default
title: meanDensityAnomalySolution
parent: WVMeanDensityAnomalyComponent
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  meanDensityAnomalySolution

return a real-valued analytical solution of the mean density anomaly mode


---

## Declaration
```matlab
 solution = meanDensityAnomalySolution(self,jMode,A,options)
```
## Parameters
+ `jMode`  integer index, (j0 >= 1 && j0 <= nModes)
+ `A`  amplitude in m.
+ `shouldAssumeConstantN`  (optional) default 1

## Returns
+ `u`  fluid velocity, u = @(x,y,z,t)
+ `v`  fluid velocity, v = @(x,y,z,t)
+ `w`  fluid velocity, w = @(x,y,z,t)
+ `eta`  isopycnal displacement, eta = @(x,y,z,t)
+ `p`  pressure, p = @(x,y,z,t)

## Discussion

  Returns function handles of the form u=@(x,y,z,t)
 
                    
