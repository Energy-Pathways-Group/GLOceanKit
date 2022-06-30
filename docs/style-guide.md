#  WaveVortexModel Style Guide

## Naming

- Class names use `CamelCase` starting with a capital letter.
- Property and method names use `camelCase` starting with a lower-case letter.
- The only exceptions are properties where the notation follows standard practice or follows a manuscript's notation. For example, the squared buoyancy frequency is written as `N2` and not `n2`.
- Properties that return a boolean should like `isPeriodic`, `shouldAntiAlias`, etc.

## Initialization

All classes can be initialized directly, or from NetCDF file and thus also written to file. 

- A class has a single initialization method (constructor), such as ` WaveVortexTransformConstantStratification(Lxyz, Nxyz, N0, options)`
- Each class has an instance method `-writeToFilePath` and `-writeToNetCDFFile`.
- Each class has a class method (or static method) to re-initialize from file, e.g., `-transformFromFilePath' and `transformFromNetCDFFile'.


