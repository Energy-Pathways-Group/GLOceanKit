---
layout: default
title: A0_QGPV_factor
parent: WVTransform
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  A0_QGPV_factor

multiplicative factor that multiplies $$A_0$$ to compute quasigeostrophic potential vorticity (QGPV).


---

## Description
Real valued transform property with dimensions $$(k,l,j)$$ and units of $$m^{-1} s^{-1}$$.

## Discussion

These coefficients multiply $$A_0$$ to give quasigeostrophic potential vorticity (QGPV).

For $$j>0$$, both the geostrophic ($$K_h>0$$) and mean-density anomaly ($$K_h=0$$) modes have coefficients computed with

$$
\textrm{QGPV}^{klj} = -\frac{g}{f_0} \left( K^2 + L_r^{-2} \right)
$$ 

where $$L_r^2 = \frac{g h}{f_0^2}$$ is the squared Rossby radius of deformation for each mode.

In the case of a rigid-lid there exists a the barotropic mode with no vortex stretching and thus, 

$$
\textrm{QGPV}^{kl0} = -\frac{g}{f_0} K^2
$$ 

for the $$j=0$$ barotropic mode.

