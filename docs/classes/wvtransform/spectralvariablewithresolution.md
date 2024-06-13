---
layout: default
title: spectralVariableWithResolution
parent: WVTransform
grand_parent: Classes
nav_order: 140
mathjax: true
---

#  spectralVariableWithResolution

create a new variable with different resolution


---

## Declaration
```matlab
 varX2 = spectralVariableWithResolution(var,Nklj)
```
## Parameters
+ `var`  a variable with dimensions [Nj Nkl]
+ `wvtX2`  a WVTransform of different size.

## Returns
+ `varX2`  matrix the size Nklj

## Discussion

  Given a variable with dimensions [Nj Nkl], this returns a new variable
  with dimensions matching wvtX2. This can be either increased or decreased
  resolution.
 
          
