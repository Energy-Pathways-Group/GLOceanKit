---
layout: default
title: masksForAllFlowConstituents
parent: WVTransform
grand_parent: Classes
nav_order: 124
mathjax: true
---

#  masksForAllFlowConstituents

Returns six 'masks' (matrices with 1s or 0s) indicating where the six


---

## Declaration
```matlab
 [IO,SGW,IGW,MDA,SG,IG] = masksForAllFlowConstituents()
```
## Returns
+ `IO`  mask for inertial oscillation solutions
+ `SGW`  mask for surface gravity wave solutions
+ `IGW`  mask for internal gravity wave solutions
+ `MDA`  mask for mean density anomaly solutions
+ `SG`  mask for surface geostrophic solutions
+ `IG`  mask for internal geostrophic solutions

## Discussion
different solution types live in the Ap, Am, A0 matrices.
 
  IO, SGW, and IGW indicate where the inertial oscillation (IO), surface
  gravity waves (SGW), and internal gravity wave (IGW) solutions live in
  the Ap and Am matrices.
 
  MDA, SG, and IG indicate where the mean density anomaly (MDA), surface
  geostrophic (SG), and interior geostrophic (IG) solutions live in the A0
  matrix.
 
  For example, if you define A = IGW .* Ap; then A will contain only the
  positive frequency internal gravity solutions.
 
                
