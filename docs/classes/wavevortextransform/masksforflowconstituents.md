---
layout: default
title: MasksForFlowConstituents
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 41
---

#  MasksForFlowConstituents

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


---

## Discussion

  Returns a sets of 'masks' (matrices with 1s or 0s) indicating where
  different solution types live in the Ap, Am, A0 matrices.
 
  Basic usage,
  [ApmMask,A0Mask] = wvm.MasksForFlowContinuents(FlowConstituents('internalGravityWave','inertialOscillation');
  will return a mask that contains 1 at the locations of the internal
  gravity waves and inertial oscillations in the Ap/Am matrices. All other
  entries will be zero.
 
  For example, if you define A = ApmMask .* Ap; then A will contain only the
  positive frequency internal gravity solutions and half the inertial solutions.
