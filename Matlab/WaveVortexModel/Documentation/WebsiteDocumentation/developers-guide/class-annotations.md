---
layout: default
title: Class annotations
parent: Developers guide
mathjax: true
---

#  Class annotations

The WaveVortexModel code uses Class Annotations to add metadata about class properties and methods. This metadata is used to support the documentation creation and writing to NetCDF file.

Matlab classes have two primary functionality groups, property and methods, which are annotated using the `CAPropertyAnnotation` and `CAMethodAnnotation`. Both of these annotations have a `detailedDescription` method that returns text from a markdown file that can be included in online documentation.

`CAPropertyAnnotation` has several useful subclasses, including `CAFunctionProperty` to annotate properties that return function handles, `CADimensionProperty` to annotate properties that should be interpreted as dimensions (such as  $$x$$, $$t$$, or $$k$$), `CANumericProperty` to annotation properties like maybe multi-dimensional arrays or scalar values, and `CAObjectProperty` which annotates properties that are array's of other objects that also inherit from `CAAnnotatedClass`.

The WaveVortexModel makes a special subclass of `CANumericProperty` called `WVVariableAnnotation` which annotated the results of `WVOperation` instances. These properties are not properties that have been declared in the class, but are computed on demand by a `WVOperation`. 

To use these annotation, a class should subclass `CAAnnotatedClass`, which provides methods for adding and removing various annotations. If the class implements `classDefinedPropertyAnnotations`,  `classRequiredPropertyNames` and initializes with options of the same name, then class instances can automatically write to a NetCDF file and be intializated from one. The `classDefinedPropertyAnnotations` function should return annotations for all relevant properties, while `classRequiredPropertyNames` should return the names of the properties that are necessary to initialize the class back to its full state. This enables `writeToFile` to work, although in practice, this it is assumed this functionally will be overridden.

Note that `-writeToFile` is an instance method, because it saves the current state of the object, while `+annotatedClassFromFile` is a static (or class) method, because it is constructing a new instance.

#  Annotations and multiple inheritance

In order to make the code more modular, the WaveVortexModel uses multiple inheritance, Matlab's approach to extending and combining class functionality. Multiple inheritance is quite tricky to work with, so I tried to establish a consistent pattern for the WaveVortexModel.

The basic class constructor pattern used for the WaveVortexModel is to have unnamed required arguments, followed by optional name-value pair arguments. For example,

```matlab
 arguments
     Lxy (1,2) double {mustBePositive}
     Nxy (1,2) double {mustBePositive}
     geomOptions.shouldAntialias (1,1) logical = true
     rotatingOptions.rotationRate (1,1) double = 7.2921E-5
     rotatingOptions.latitude (1,1) double = 33
     rotatingOptions.g (1,1) double = 9.81
     options.h (1,1) double = 0.8
 end
```

requires the size and number of grid points, but everything else can be set by a default value. This works when creating a new model, however, to recover an instance with the exact same *state* requires all properties, including optional properties, to be specified. Hence, this class would need to list `Lx`, `Ly`, `Nx`, `Ny`, `shouldAntialias`, etc. as `+classRequiredPropertyNames` in order to save and then restore the state from those saved values. Additionally, because the class constructor follows a specialized pattern, the class needs to define a static method `+waveVortexTransformFromFile()` with a custom implementation.

Where this gets complicated is with multiple inheritance. A class inheriting annotations from multiple classes is one thing, but then it also has to initialize each of the super classes correctly. The rules are as follows:

All classes have the name structure `WV[ClassName]`, and in fact use `WV[RootClassName][SubclassName]` to extend functionality.

## Root class annotations

 So each root class and subclass should implement,

- `+propertyAnnotationsFor[RootClassName]` to return all property annotations for the class, and
- `+namesOfRequiredPropertiesFor[RootClassName]` to return the names of all required properties.

And thus implementing,

```matlab
function propertyAnnotations = classDefinedPropertyAnnotations()
    propertyAnnotations = WV[RootClassName].propertyAnnotationsFor[RootClassName]();
end

function vars = classRequiredPropertyNames()
    vars = WV[RootClassName].namesOfRequiredPropertiesForStratification[RootClassName]();
end
```
 
will satisfy the requirements of the `CAAnnotatedClass` with the correct annotations.

To support intialization, the root class must read the properties from file and call the class constructor. All class therefore have the following methods,

- `+requiredPropertiesFor[RootClassName]FromGroup` which returns properties that match the signature of the class constructor, e.g., `[Lxy,Nxy,options]` for the example above.
- `+[rootClassName]FromGroup` which calls the above, and returns a newly constructed instance, and finally
- `+[rootClassName]FromFile` which calls the above, and returns a newly constructed instance.

Implementation of `+[rootClassName]FromGroup` should use `CAAnnotatedClass.canInitializeDirectlyFromGroup` or `CAAnnotatedClass.throwErrorIfMissingProperties` to confirm that the required properties are present, followed by `CAAnnotatedClass.propertyValuesFromGroup` to fetch the properties. 

These are all fairly trivial to implement using the supported functions, but they must be implemented.

## Subclass annotations

A subclass might

1. gain new annotated properties and new required properties,
2. lose required properties (because the subclass makes assumptions/forces values),
3. inherit from two or more annotated root classes.



```matlab
function propertyAnnotations = classDefinedPropertyAnnotations()
    propertyAnnotations = WV[RootClassName][SubclassName].propertyAnnotationsFor[RootClassName]();
end

function vars = classRequiredPropertyNames()
    vars = WV[RootClassName][SubclassName].namesOfRequiredPropertiesForStratification[RootClassName]();
end
```
