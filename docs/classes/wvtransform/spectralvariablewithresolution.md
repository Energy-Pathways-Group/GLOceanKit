---
layout: default
title: spectralVariableWithResolution
parent: WVTransform
grand_parent: Classes
nav_order: 177
mathjax: true
---

#  spectralVariableWithResolution

create a new variable with increased resolution


---

## Declaration
```matlab
 varX2 = spectralVariableWithResolution(var,Nklj)
```
## Parameters
+ `var`  a variable with dimensions [Nk Nl Nj]
+ `Nklj`  vector of size [1 3] with new dimensions, [NkX2 NlX2 NjX2]

## Returns
+ `varX2`  matrix the size Nklj

## Discussion

  Given a variable with dimensions [Nk Nl Nj], this returns a new variable
  with the dimension Nklj.
 
          
