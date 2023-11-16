---
layout: default
title: A0_PE_factor
parent: WVTransform
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  A0_PE_factor

multiplicative factor that multiplies $$A_0^2$$ to compute potential energy.


---

## Description
Real valued transform property with dimensions $$(k,l,j)$$ and units of $$m s^{-2}$$.

## Discussion

These coefficients multiply $$A_0^2$$ to give a horizontally-averaged depth-integrated potential energy for the geostrophic solutions.

Assuming hydrostatic modes, this is

$$
PE = \frac{g}{2}
$$ 

for the $$j>0$$ modes and

$$
PE = 0
$$ 

for the $$j=0$$ barotropic mode.

