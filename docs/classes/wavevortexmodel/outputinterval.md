---
layout: default
title: outputInterval
parent: WaveVortexModel
grand_parent: Classes
nav_order: 10
---

#  outputInterval

Model output interval (seconds)


---

## Discussion
This property is optionally set when calling setupIntegrator. If
  set, it will allow you to call -integrateToNextOutputTime and, if
  a NetCDF file is set for output, it will set the interval at
  which time steps are written to file.