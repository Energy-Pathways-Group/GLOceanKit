---
layout: default
title: A0_TE_factor
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  A0_TE_factor

multiplicative factor that multiplies $$A_0$$ to compute total energy.


---

## Description
Real valued transform property with dimensions $$(k,l,j)$$ and units of $$m s^{-2}$$.

## Discussion

These coefficients multiply $$A_0$$ to give a horizontally-averaged depth-integrated total energy for the geostrophic solutions.

Assuming hydrostatic modes, this is

$$
E = \frac{g}{2} \left( \frac{gh}{f_0^2} K^2 + 1 \right)
$$ 

and

$$
E = \frac{g}{2} \left( \frac{g L_z}{f_0^2} K^2 + 1 \right)
$$ 

for the barotropic mode.
