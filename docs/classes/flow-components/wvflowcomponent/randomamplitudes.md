---
layout: default
title: randomAmplitudes
parent: WVFlowComponent
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  randomAmplitudes

returns random amplitude for a valid flow state


---

## Declaration
```matlab
 Ap,Am,A0] = randomAmplitudes()
```
## Returns
+ `Ap`  matrix of size [Nj Nkl]
+ `Am`  matrix of size [Nj Nkl]
+ `A0`  matrix of size [Nj Nkl]

## Discussion

  Returns Ap, Am, A0 matrices initialized with random amplitude
  for this flow component. These resulting matrices will have
  the correct symmetries for a valid flow state. 
 
          
