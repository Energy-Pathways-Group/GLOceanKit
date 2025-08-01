---
layout: default
title: Wavenumbers, modes, and indices
parent: Users guide
mathjax: true
nav_order: 6
---

#  Wavenumbers, modes, and indices

Working with the WaveVortexModel requires working with wavenumbers (such as $$k$$ and $$l$$), modes (such as geostrophic modes, or mode number $$j$$), and indices into matrices. There are many best practices and tools for working with wavenumbers, mode, and indices.

## Definitions

A **wavenumber** is a [spatial frequency](https://en.wikipedia.org/wiki/Wavenumber) with units of radians per meter. All WVTransform classes have horizontal wavenumbers $$k$$ and $$l$$, and the constant stratification class also has a vertical wavenumber $$m$$.

A **mode** is a unitless integer quantity. The model refers to `kMode` and `lMode`, and `j` all as unitless modes. Additionally, a geostrophic solution is also a mode (as are the other linear solutions).

An **index** is a location in memory of a particular variable. Matlab uses [subscript indices, linear indices, and logical indexing](https://www.mathworks.com/company/technical-articles/matrix-indexing-in-matlab.html) for different ways of accessing memory. The prefered method of indexing is local indexing.

## More about modes

Each degree-of-freedom in the wave-vortex model is a mode, which, in this case, is specifically a linear solution to the equations of motion. Each degree-of-freedom is identified with a unique set of mode numbers, either `(kMode,lMode)` or `(kMode,lMode,jMode)`, depending on whether the flow is two-dimensional or three-dimensional. There are two aspects of this identification which are important to call out:
1. each set of mode numbers identifies a unique solution, and
2. the choice of mode number for each solution is a matter of convention.

As a consequence of these two aspects, you, the user, cannot know *a priori* what a valid set of mode numbers is unless the notation is established. The terminology used here is that these mode numbers are *primary* mode numbers. Furthermore, there may be other reasonable choices of mode number which are standard and perfectly valid. The terminology used here is that these are *conjugate* mode numbers. If a conjugate mode number exists, then this implies that there is some mapping from the conjugate mode to the primary solution. For these reasons, the all `WVTransform` subclass respond to the follow APIs,

- `bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode)`
- `bool = isValidConjugateModeNumber(self,kMode,lMode,jMode)`
- `bool = isValidModeNumber(self,kMode,lMode,jMode)`

The function `isValidPrimaryModeNumber` will tell you whether or not that particular set of mode numbers is a valid, unique solution that is part of the basis set being used. The function `isValidConjugateModeNumber` will similarly tell you whether that mode number corresponds to a valid solution, but one which is not primary.

As a concrete example, consider the the modal solution to the periodic harmonic equation, $$u = A \sin\left( n \frac{\pi x}{L} \right)$$. A resonable convention would be to establish that $$n$$ must be an integer such that $$ 1 \leq n \leq 7$$ as the 7 primary mode numbers. Alternatively, it would be perfectly reasonable to choose $$-7 \leq n \leq -1$$ as the mode numbers and thus as a matter of convention we could identify these as conjugate mode numbers. Both are valid sets of mode numbers, and the mapping that exists between these different set of mode numbers requires we change the sign on the amplitude of the solution.

The situation gets more complicated as the solutions get more complicated. For doubly periodic geometry we now make a choice as to which set of wavenumber in which dimension gets designated as conjugate, and thus `WVGeometryDoublyPeriodic` has a `conjugateDimension` property which gets set upon initialization. Use of the above APIs let this choice be invisible to the user, even though it strongly affects the numerics.

The modes of the `WVTransform` classes are not solutions to the harmonic equation, but are instead solutions to various linearized equations of motion of geophysical flows. These solutions include the usual geostrophic and internal gravity solutions. The wave-vortex model identifies these solution types using the `WVFlowComponent` subclass `WVPrimaryFlowComponent`. The flow components also respond to the same three APIs listed above.
