---
layout: default
title: WaveVortexModel
parent: Classes
---
#  WaveVortexModel

The WaveVortexModel is responsible for time-stepping (integrating) the ocean state forward in time using a WaveVortexTransform.

## Overview



## Topics
+ [Initialization](#initialization)
  + [`WaveVortexModel`](#dimensionwithname)
+ [Particles](#ocean-state)
  + [`AddParticles`](#attributes)
  + [`ParticlePositions`](#attributes)
+ [Tracer](#ocean-state)
  + [`AddTracer`](#attributes)
  + [`Tracer`](#attributes)
+ [Reduced Interaction Models](#ocean-state)
  + Limiting interactions
    + [`allowNonlinearInteractionsWithModes`](#attributes)
    + [`allowNonlinearInteractionsWithConstituents`](#attributes)
    + [`disallowNonlinearInteractionsWithConstituents`](#attributes)
  + Limiting energy fluxes
    + [`unfreezeEnergyOfConstituents`](#attributes)
    + [`freezeEnergyOfConstituents`](#attributes)
  + Anti-aliasing
    + [`IsAntiAliased`](#attributes)
    + [`disallowNonlinearInteractionsWithAliasedModes`](#attributes)
    + [`freezeEnergyOfAliasedModes`](#attributes)
    + [`clearEnergyFromAliasedModes`](#attributes)
  + Forcing
    + [`addForcingWaveModes`](#attributes)
+ [Integration](#ocean-state)
  + [`SetupIntegrator`](#attributes)
  + [`IntegrateToTime`](#attributes)
  + [`integrateOneTimeStep`](#attributes)
  + [`integrateToNextOutputTime`](#attributes)
+ [Writing To File](#ocean-state)
  + [`CreateNetCDFFileForModelOutput`](#attributes)
