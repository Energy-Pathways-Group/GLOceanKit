---
layout: default
title: masksForAllFlowConstituents
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 120
mathjax: true
---

#  masksForAllFlowConstituents

Returns six 'masks' (matrices with 1s or 0s) indicating where the six


---

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
 
  
