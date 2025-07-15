---
layout: default
title: Reading and writing to file, advanced topics
parent: Users guide
mathjax: true
nav_order: 3
has_toc: true
---

#  Reading and writing to file, advanced topics

The `WVModel` supports a number of advanced features for writing to file, including
- multiple output files
- groups within files with different output intervals
- groups that start and stop at different time points

With these features, it is possible to setup up numerical simulations that "observe" the ocean in different ways, at different intervals. For example, you might place moorings in your simulation that sample the velocity field at a high frequency. You might setup a series of drifter experiments, that deploy and retrieve drifters every five days.

There are three classes that work together to write to file,
- `WVModelObservingSystem` - classes that describe different ways of observing the fluid, may or may not require integration
- `WVModelOutputFile` - a representation of a file to be written to disk; has one more output groups
- `WVModelOutputGroup` - a netcdf group that writes to file at certain output times; has one or more observing systems

Any `WVModelObservingSystem` instances that require integration, will also be held onto the model, which will integrate the observign system as necessary.


## WVModelObservingSystem

Observing systems include, e.g., the wave-vortex coefficients, Eulerian fields, Lagrangian particles, tracers, mooring, and satellite along-track data.

Subclasses of `WVModelObservingSystem` describe whether or not the observing system needs to be integrated (fluxed) in time by the model, and how to write to a group (an instance of `WVModelOutputGroup`) if desired. Some observing seems are only fluxed (e.g., the `WVCoefficients` system is used to integrate the wave-vortex coefficients) and do not write to file, while other observing system are not fluxed (e.g., `WVMooring`) but do write to file. Some of the observing systems require a specific subclasses of `WVModelOutputGroup` in order to write to file and set `requiresCustomOutputTimes` to `true`. For example, the `WVAlongTrackObservingSystem` can only write to the `WVModelOutputGroupAlongTrack` because it requires specific output times that are tied to the details of the satellite oribits.

## WVModelOutputFile

After you create a `WVModel` instance, you may add one more `WVModelOutputFile`. Output files are conceptually very simple as they simply represent files on disk. Importantly, they hold onto one or more output groups, instances of `WVModelOutputGroup`, and internally they orchestrate pausing the model and writing to groups.

## WVModelOutputGroup

The most useful `WVModelOutputGroup` subclass is the `WVModelOutputGroupEvenlySpaced`, which outputs at a fixed interval, but with user specified initial and final times. 
