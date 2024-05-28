---
layout: default
title: throwErrorIfDensityViolation
parent: WVTransformBoussinesq
grand_parent: Classes
nav_order: 30
mathjax: true
---

#  throwErrorIfDensityViolation

checks if the proposed coefficients are a valid adiabatic re-arrangement of the base state


---

## Discussion

  Given some proposed new set of values for A0, Ap, Am, will
  the fluid state violate our density condition? If yes, then
  throw an error and tell the user about it.
