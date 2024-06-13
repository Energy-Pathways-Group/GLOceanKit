---
layout: default
title: operationForDynamicalVariable
parent: WVTransform
grand_parent: Classes
nav_order: 123
mathjax: true
---

#  operationForDynamicalVariable

This function is designed with the following goals:


---

## Discussion
1) the definition of each variable is made *once*, with the intention
    of minimizing mistakes and making it easier to read. The challenge with
    the goals is then that masks have to be applied generically.
    2) Masking an entire component (Apm or A0) to zero should avoid any
    extra computation. This is important because we still want the QG flow
    to retain its speed.
 
  The solution I came up with was that each variables must specify the
  recipe for computing itself as three function handles that tell it how to
  construct the three components Ap,Am,A0. The function then only applies
  the function handle when necessary.
