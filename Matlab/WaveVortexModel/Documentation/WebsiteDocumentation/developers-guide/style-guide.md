---
layout: default
title: WaveVortexModel style guide
parent: Developers guide
mathjax: true
---

#  WaveVortexModel style guide

## Naming

- All classes within the WaveVortexModel use a prefix of `WV` to establish a namespace.
- Any user defined subclasses should *not* use the `WV` prefix.
- Class names use [`UpperCamelCase`](https://en.wikipedia.org/wiki/Camel_case) starting with a capital letter.
- Property and method names use [`lowerCamelCase`](https://en.wikipedia.org/wiki/Camel_case) starting with a lower-case letter.
- The only exceptions are properties where the notation follows standard practice or follows a manuscript's notation. For example, the squared buoyancy frequency is written as `N2` and not `n2`.
- Properties that return a boolean should read like an assertion, e.g., `isPeriodic` or `shouldAntiAlias`.

## Initialization

All classes can be initialized directly, or from NetCDF file and thus also written to file. 

- A class has a single initialization method (constructor), such as ` WaveVortexTransformConstantStratification(Lxyz, Nxyz, N0, options)`
- Each class has an instance method `-writeToFilePath` and `-writeToNetCDFFile`.
- Each class has a class method (or static method) to re-initialize from file, e.g., `-transformFromFilePath' and `transformFromNetCDFFile'.


