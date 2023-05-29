---
layout: default
title: WV transform & model
nav_order: 4
has_children: true
permalink: /classes-transform-and-model
---

#  Classes

The `WaveVortexModel` consists of several classes and subclasses.

- [WVTransform](/classes/wvtransform/)
- [WVModel](/classes/wvmodel/)
- [NetCDFFile](/classes/netcdffile/)

## WVTransform

The WVTransform subclasses encapsulate data representing the *state* of the ocean at a given instant in time (e.g., u, v, w, and rho). What makes the WaveVortexTransform subclasses special is that the state of the ocean is represented as energetically independent waves and geostrophic motions (vortices). These classes can be queried for other information, e.g., Ertel PV, relative vorticity, etc.

- [WVTransformConstantStratification]()
- [WVTransformHydrostatic]()

## WVModel

The WaveVortexModel uses the WaveVortexTransform to integrate (time-step) the non-linear equations of motion forward in time. The model adds robust support for particle advection, tracer advection, as well as reduced interaction models.

## NetCDFFile

Tools for reading and writing the model to file, including support for tracers, floats and model restarts.

- [NetCDFFile](/classes/netcdffile/)
