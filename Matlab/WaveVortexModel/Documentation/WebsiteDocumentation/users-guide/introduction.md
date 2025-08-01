---
layout: default
title: Introduction
parent: Users guide
mathjax: true
nav_order: 1
---

#  Introduction

There are two high-level components to the wave-vortex model, the `WVTransform` and the `WVModel`. The `WV` prefix is short for wave-vortex and establishes a namespace for all classes in the model.

## WVTransform

The `WVTransform` subclasses encapsulate data representing the *state* of the ocean at a given instant in time (e.g., $$u$$, $$v$$, $$w$$, and $$\rho$$). The `WVTransform` can be customized to tell you anything you want about the state of the ocean at a point in time.

## WVModel

A `WVModel` instance uses a `WVTransform` instance to integrate (time-step) the non-linear equations of motion forwad in time. The `WVModel` class adds support for features like particle advection and tracer advection.


