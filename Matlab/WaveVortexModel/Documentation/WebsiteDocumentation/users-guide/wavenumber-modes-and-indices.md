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



