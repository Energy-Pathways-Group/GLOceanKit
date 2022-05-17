---
layout: default
title: Classes
nav_order: 2
has_children: true
permalink: docs/classes
---

#  Classes

The `WaveVortexModel` consists of several classes and subclasses.

- [WaveVortexModel]()
- [WaveVortexModelIntegrationTools]()
- [WaveVortexModelNetCDFFile]()

## WaveVortexModel

There are two usable subclasses at the moment,

- [WaveVortexModelConstantStratification]()
- [WaveVortexModelHydrostatic]()

## WaveVortexModelIntegrationTools

Tools for integrating (time-stepping) the model forward.

## WaveVortexModelNetCDFFile

Tools for reading and writing the model to file, including support for tracers, floats and model restarts.

- [NetCDFFile](./netcdffile.html)