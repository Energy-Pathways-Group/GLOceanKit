---
layout: default
title: A0_TZ_factor
parent: WVTransform
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  A0_TZ_factor

multiplicative factor that multiplies $$A_0^2$$ to compute quasigeostrophic enstrophy.


---

## Description
Real valued transform property with dimensions $$(k,l,j)$$ and units of $$m^{-1} s^{-2}$$.

## Discussion

These coefficients multiply $$A_0^2$$ to give a horizontally-averaged depth-integrated total quasigeostrophic enstrophy.

For $$j>0$$, both the geostrophic ($$K_h>0$$) and mean-density anomaly ($$K_h=0$$) modes have coefficients computed with

$$
\textrm{TZ}^{klj} = \frac{g}{2} \left( K^2 + L_r^{-2} \right)^2 L_r^2 
$$ 

where $$L_r^2 = \frac{g h}{f_0^2}$$ is the squared Rossby radius of deformation for each mode. In the case of a rigid-lid there exists a the barotropic mode with no vortex stretching and thus, 

$$
\textrm{TZ}^{kl0} = \frac{g}{2} K^4 \frac{g L_z}{f_0^2}
$$ 

for the $$j=0$$ barotropic mode.

