---
layout: default
title: A0_HKE_factor
parent: WVTransform
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  A0_HKE_factor

multiplicative factor that multiplies $$A_0^2$$ to compute horizontal kinetic energy.


---

## Description
Real valued transform property with dimensions $$(k,l,j)$$ and units of $$m s^{-2}$$.

## Discussion

These coefficients multiply $$A_0^2$$ to give a horizontally-averaged depth-integrated horizontal kinetic energy for the geostrophic solutions.

Assuming hydrostatic modes, this is

$$
HKE = \frac{g}{2} \left( \frac{gh}{f_0^2} K^2 \right)
$$ 

for the $$j>0$$ modes and

$$
HKE = \frac{g}{2} \left( \frac{g L_z}{f_0^2} K^2\right)
$$ 

for the $$j=0$$ barotropic mode.

