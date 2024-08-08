---
layout: default
title: Annotations
parent: Developers guide
mathjax: true
---

#  Annotations

The methods and properties of the `WVTransform' are divided into 5 categories, depending on their functionality. These properties or methods are then annotated appropriately.

- *Dimension* are independent coordinate axes, such as $$x$$, $$t$$, or $$k$$ and are annotated with `WVDimensionAnnotation`.
- *Properties* are fixed attributes of the `WVTransform` that are never time-dependent. The `WVPropertyAnnotation` class is used to annotate properties.
- *Variables* depend on the current state of flow and always depend on the wave-vortex coefficients $$(A_+, A_-, A_0)$$. Variables are annotated with `WVVariableAnnotation`.

Each dimension, property, or variable can be individually access from a `WVTransform` instance, e.g., `wvt.u` will return the `u` variable. However, there are also methods which can return multiple dimensions, properties, or variables, e.g., `[u,v,w] = wvt.velocityField` is a method that returns three variables.

We need to distinguish between variables that have direct method (or even property) calls, and those which need to be computed, and thus might be stored in the cache.

- *Methods* are everything else, and use `WVAnnotation` to annotate their functionality.

If you are building a subclass the `WVTransform`, it is important to know that there are two types of variables: those that are directly computed by the class, and those that result from a `WVOperation`.

