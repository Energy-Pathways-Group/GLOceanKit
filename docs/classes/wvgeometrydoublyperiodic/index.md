---
layout: default
title: WVGeometryDoublyPeriodic
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 7
---

#  WVGeometryDoublyPeriodic

A domain periodic in both x and y.


---

## Overview
 
  The WVGeometryDoublyPeriodic encapsulates the transformations and
  logic necessary for computing Fourier transforms and spectral
  derivatives in a doubly periodic domain.
 
  The class primarily converts between two different data structures
  for representing functions of (k,l): what we call the "DFT" grid, and
  the "WV" grid. The DFT grid is appropriate for many FFT algorithms,
  while the WV grid is ideal for the vertical mode matrix
  multiplications in the WaveVortexModel, as it contains no redundant
  coefficients.
 
  The DFT layout is the the matrix structure used by many modern FFT
  algorithms. For real-valued functions, this format contains twice as
  much information as necessary, due to Hermitian conjugacy.
  Additionally, we often want to restrict ourselves to wavenumbers that
  do not alias with quadratic multiplication (the two-thirds rule). In
  two dimensions, this means that 4/9ths of the available wavenumbers
  are aliased. The WV layout thus does not include either the Hermitian
  conjugates nor the aliased modes.
 
  The basic usage of the indices is as follows:
  assume wvMatrix and dftMatrix are shaped as
    size(wvMatrix) = [Nkl_wv 1]
    size(dftMatrix) = [Nk_dft Nl_dft] (equivalently [Nx Ny]
  then to transform data from the DFT matrix to the WV matrix,
    wvMatrix = dftMatrix(dftPrimaryIndices);
  and the reverse is
    dftMatrix(dftPrimaryIndices) = wvMatrix;
    dftMatrix(dftConjugateIndices) = conj(wvMatrix(wvConjugateIndex));


## Topics
+ Other
  + [`Lx`](/classes/wvgeometrydoublyperiodic/lx.html) 
  + [`Ly`](/classes/wvgeometrydoublyperiodic/ly.html) 
  + [`Nk_dft`](/classes/wvgeometrydoublyperiodic/nk_dft.html) 
  + [`Nkl_wv`](/classes/wvgeometrydoublyperiodic/nkl_wv.html) 
  + [`Nl_dft`](/classes/wvgeometrydoublyperiodic/nl_dft.html) 
  + [`Nx`](/classes/wvgeometrydoublyperiodic/nx.html) 
  + [`Ny`](/classes/wvgeometrydoublyperiodic/ny.html) 
  + [`conjugateDimension`](/classes/wvgeometrydoublyperiodic/conjugatedimension.html) 
  + [`dftConjugateIndices`](/classes/wvgeometrydoublyperiodic/dftconjugateindices.html) 
  + [`dftPrimaryIndices`](/classes/wvgeometrydoublyperiodic/dftprimaryindices.html) 
  + [`kMode_dft`](/classes/wvgeometrydoublyperiodic/kmode_dft.html) 
  + [`kMode_wv`](/classes/wvgeometrydoublyperiodic/kmode_wv.html) 
  + [`k_dft`](/classes/wvgeometrydoublyperiodic/k_dft.html) 
  + [`k_wv`](/classes/wvgeometrydoublyperiodic/k_wv.html) 
  + [`lMode_dft`](/classes/wvgeometrydoublyperiodic/lmode_dft.html) 
  + [`lMode_wv`](/classes/wvgeometrydoublyperiodic/lmode_wv.html) 
  + [`l_dft`](/classes/wvgeometrydoublyperiodic/l_dft.html) 
  + [`l_wv`](/classes/wvgeometrydoublyperiodic/l_wv.html) 
  + [`shouldAntialias`](/classes/wvgeometrydoublyperiodic/shouldantialias.html) 
  + [`shouldExcludeNyquist`](/classes/wvgeometrydoublyperiodic/shouldexcludenyquist.html) 
  + [`shouldExludeConjugates`](/classes/wvgeometrydoublyperiodic/shouldexludeconjugates.html) 
  + [`transformFromSpatialDomain`](/classes/wvgeometrydoublyperiodic/transformfromspatialdomain.html) 
  + [`transformToSpatialDomain`](/classes/wvgeometrydoublyperiodic/transformtospatialdomain.html) 
  + [`x`](/classes/wvgeometrydoublyperiodic/x.html) 
  + [`y`](/classes/wvgeometrydoublyperiodic/y.html) 
+ Initialization
  + [`WVGeometryDoublyPeriodic`](/classes/wvgeometrydoublyperiodic/wvgeometrydoublyperiodic.html) create a geometry for a  doubly periodic domain
+ Utility function
  + [`degreesOfFreedomForComplexMatrix`](/classes/wvgeometrydoublyperiodic/degreesoffreedomforcomplexmatrix.html) a matrix with the number of degrees-of-freedom at each entry
  + [`degreesOfFreedomForRealMatrix`](/classes/wvgeometrydoublyperiodic/degreesoffreedomforrealmatrix.html) a matrix with the number of degrees-of-freedom at each entry
  + [`indicesOfFourierConjugates`](/classes/wvgeometrydoublyperiodic/indicesoffourierconjugates.html) a matrix of linear indices of the conjugate
  + [`isHermitian`](/classes/wvgeometrydoublyperiodic/ishermitian.html) Check if the matrix is Hermitian. Report errors.
  + [`setConjugateToUnity`](/classes/wvgeometrydoublyperiodic/setconjugatetounity.html) set the conjugate of the wavenumber (iK,iL) to 1
+ Operations
  + Differentiation
    + [`diffX`](/classes/wvgeometrydoublyperiodic/diffx.html) differentiate a spatial variable in the x-direction
    + [`diffY`](/classes/wvgeometrydoublyperiodic/diffy.html) differentiate a spatial variable in the y-direction
+ Index Gymnastics
  + [`indicesFromDFTGridToWVGrid`](/classes/wvgeometrydoublyperiodic/indicesfromdftgridtowvgrid.html) indices to convert from DFT to WV grid
  + [`indicesFromWVGridToDFTGrid`](/classes/wvgeometrydoublyperiodic/indicesfromwvgridtodftgrid.html) indices to convert from WV to DFT grid
  + [`isValidWVModeNumber`](/classes/wvgeometrydoublyperiodic/isvalidwvmodenumber.html) return a boolean indicating whether (k,l) is a valid WV mode number
  + [`linearWVIndexFromModeNumber`](/classes/wvgeometrydoublyperiodic/linearwvindexfrommodenumber.html) return the linear index into k_wv and l_wv from a mode number
  + [`modeNumberFromWVIndex`](/classes/wvgeometrydoublyperiodic/modenumberfromwvindex.html) return mode number from a linear index into a WV matrix
  + [`transformFromDFTGridToWVGrid`](/classes/wvgeometrydoublyperiodic/transformfromdftgridtowvgrid.html) convert from DFT to WV grid
  + [`transformFromWVGridToDFTGrid`](/classes/wvgeometrydoublyperiodic/transformfromwvgridtodftgrid.html) convert from a WV to DFT grid
+ Masks
  + [`maskForAliasedModes`](/classes/wvgeometrydoublyperiodic/maskforaliasedmodes.html) returns a mask with locations of modes that will alias with a quadratic multiplication.
  + [`maskForConjugateFourierCoefficients`](/classes/wvgeometrydoublyperiodic/maskforconjugatefouriercoefficients.html) a mask indicate the components that are redundant conjugates
  + [`maskForNyquistModes`](/classes/wvgeometrydoublyperiodic/maskfornyquistmodes.html) returns a mask with locations of modes that are not fully resolved


---