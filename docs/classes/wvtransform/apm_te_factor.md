---
layout: default
title: Apm_TE_factor
parent: WVTransform
grand_parent: Classes
nav_order: 18
mathjax: true
---

#  Apm_TE_factor

multiplicative factor that multiplies $$A_\pm$$ to compute total energy.


---

## Description
Real valued transform property with dimensions $$(k,l,j)$$ and units of $$m$$.

## Discussion

These coefficients multiply $$A_+$$ and $$A_-$$ to give a horizontally-averaged depth-integrated total energy for those wave solutions.

For the internal waves this is simply,

$$
E = A_\pm^2 h
$$ 

and

$$
E = A_\pm^2 L_z
$$ 

for the barotropic mode at $$j=0$$.

