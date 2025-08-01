---
layout: default
title: Design of the primary flow components
parent: Developers guide
mathjax: true
nav_order: 2
---

#  Design of the primary flow components

Each primary flow component has exactly two classes to support its logic: a subclass of `WVPrimaryFlowComponent`, e.g. `WVInternalGravityWaveComponent` and a `Methods` class, e.g., `WVInternalGravityWaveMethods`. 

Any subclass that wants to support internal gravity waves must be a subclass of `WVInternalGravityWaveMethods` and then after `WVStratifiedFlow` and `WVTransform` have had their constructors called, the instance method `-initializeInternalGravityWaveComponent` can be called. This will add the primary flow component, and initialize the various variables need to support those components.

 The motivation for this design choice was to keep the logic surrounding particular flow components in one place (okay, fine; so it ended up being two places). 
