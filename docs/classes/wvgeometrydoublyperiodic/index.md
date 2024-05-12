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
  ```matlab
    size(wvMatrix) == [Nkl_wv 1];
    size(dftMatrix) == [Nk_dft Nl_dft]; % (equivalently [Nx Ny]
  ```
  then to transform data from the DFT matrix to the WV matrix,
  ```matlab
    wvMatrix = dftMatrix(dftPrimaryIndices);
  ```
  and the reverse is
  ```matlab
    dftMatrix(dftPrimaryIndices) = wvMatrix;
    dftMatrix(dftConjugateIndices) = conj(wvMatrix(wvConjugateIndex));
  ```
 
                        


## Topics
+ Initialization
  + [`WVGeometryDoublyPeriodic`](/classes/wvgeometrydoublyperiodic/wvgeometrydoublyperiodic.html) create a geometry for a  doubly periodic domain
+ Domain attributes
  + Spatial grid
    + [`Lx`](/classes/wvgeometrydoublyperiodic/lx.html) length of the x-dimension
    + [`Ly`](/classes/wvgeometrydoublyperiodic/ly.html) length of the y-dimension
    + [`Nx`](/classes/wvgeometrydoublyperiodic/nx.html) number of grid points in the x-dimension
    + [`Ny`](/classes/wvgeometrydoublyperiodic/ny.html) number of grid points in the y-dimension
    + [`x`](/classes/wvgeometrydoublyperiodic/x.html) dimension
    + [`y`](/classes/wvgeometrydoublyperiodic/y.html) dimension
  + DFT grid
    + [`Nk_dft`](/classes/wvgeometrydoublyperiodic/nk_dft.html) length of the k-wavenumber dimension on the DFT grid
    + [`Nl_dft`](/classes/wvgeometrydoublyperiodic/nl_dft.html) length of the l-wavenumber dimension on the DFT grid
    + [`conjugateDimension`](/classes/wvgeometrydoublyperiodic/conjugatedimension.html) assumed conjugate dimension
    + [`kMode_dft`](/classes/wvgeometrydoublyperiodic/kmode_dft.html) k mode-number on the DFT grid
    + [`k_dft`](/classes/wvgeometrydoublyperiodic/k_dft.html) k wavenumber dimension on the DFT grid
    + [`lMode_dft`](/classes/wvgeometrydoublyperiodic/lmode_dft.html) l mode-number on the DFT grid
    + [`l_dft`](/classes/wvgeometrydoublyperiodic/l_dft.html) l wavenumber dimension on the DFT grid
  + WV grid
    + [`Nkl_wv`](/classes/wvgeometrydoublyperiodic/nkl_wv.html) length of the combined kl-wavenumber dimension on the WV grid
    + [`dftConjugateIndices`](/classes/wvgeometrydoublyperiodic/dftconjugateindices.html) index into the DFT grid of the conjugate of each WV mode
    + [`dftPrimaryIndices`](/classes/wvgeometrydoublyperiodic/dftprimaryindices.html) index into the DFT grid of each WV mode
    + [`kMode_wv`](/classes/wvgeometrydoublyperiodic/kmode_wv.html) k mode number on the WV grid
    + [`k_wv`](/classes/wvgeometrydoublyperiodic/k_wv.html) k-wavenumber dimension on the WV grid
    + [`lMode_wv`](/classes/wvgeometrydoublyperiodic/lmode_wv.html) l mode number on the WV grid
    + [`l_wv`](/classes/wvgeometrydoublyperiodic/l_wv.html) l-wavenumber dimension on the WV grid
    + [`shouldAntialias`](/classes/wvgeometrydoublyperiodic/shouldantialias.html) whether the WV grid includes quadratically aliased wavenumbers
    + [`shouldExcludeNyquist`](/classes/wvgeometrydoublyperiodic/shouldexcludenyquist.html) whether the WV grid includes Nyquist wavenumbers
    + [`shouldExludeConjugates`](/classes/wvgeometrydoublyperiodic/shouldexludeconjugates.html) whether the WV grid includes wavenumbers that are Hermitian conjugates
+ Operations
  + Grid transformation
    + [`transformFromDFTGridToWVGrid`](/classes/wvgeometrydoublyperiodic/transformfromdftgridtowvgrid.html) convert from DFT to WV grid
    + [`transformFromWVGridToDFTGrid`](/classes/wvgeometrydoublyperiodic/transformfromwvgridtodftgrid.html) convert from a WV to DFT grid
  + Fourier transformation
    + [`transformFromSpatialDomain`](/classes/wvgeometrydoublyperiodic/transformfromspatialdomain.html) transform from $$(x,y,z)$$ to $$(k,l,z)$$ on the DFT grid
    + [`transformToSpatialDomain`](/classes/wvgeometrydoublyperiodic/transformtospatialdomain.html) transform from $$(k,l,z)$$ on the DFT grid to $$(x,y,z)$$
  + Differentiation
    + [`diffX`](/classes/wvgeometrydoublyperiodic/diffx.html) differentiate a spatial variable in the x-direction
    + [`diffY`](/classes/wvgeometrydoublyperiodic/diffy.html) differentiate a spatial variable in the y-direction
+ Index gymnastics
  + [`indicesFromDFTGridToWVGrid`](/classes/wvgeometrydoublyperiodic/indicesfromdftgridtowvgrid.html) indices to convert from DFT to WV grid
  + [`indicesFromWVGridToDFTGrid`](/classes/wvgeometrydoublyperiodic/indicesfromwvgridtodftgrid.html) indices to convert from WV to DFT grid
  + [`isValidConjugateWVModeNumber`](/classes/wvgeometrydoublyperiodic/isvalidconjugatewvmodenumber.html) return a boolean indicating whether (k,l) is a valid conjugate WV mode number
  + [`isValidPrimaryWVModeNumber`](/classes/wvgeometrydoublyperiodic/isvalidprimarywvmodenumber.html) return a boolean indicating whether (k,l) is a valid primary (non-conjugate) WV mode number
  + [`isValidWVModeNumber`](/classes/wvgeometrydoublyperiodic/isvalidwvmodenumber.html) return a boolean indicating whether (k,l) is a valid WV mode number
  + [`modeNumberFromWVIndex`](/classes/wvgeometrydoublyperiodic/modenumberfromwvindex.html) return mode number from a linear index into a WV matrix
  + [`primaryModeNumberFromWVModeNumber`](/classes/wvgeometrydoublyperiodic/primarymodenumberfromwvmodenumber.html) takes any valid WV mode number and returns the primary mode number
  + [`wvIndexFromModeNumber`](/classes/wvgeometrydoublyperiodic/wvindexfrommodenumber.html) return the linear index into k_wv and l_wv from a mode number
+ Masks
  + [`maskForAliasedModes`](/classes/wvgeometrydoublyperiodic/maskforaliasedmodes.html) returns a mask with locations of modes that will alias with a quadratic multiplication.
  + [`maskForConjugateFourierCoefficients`](/classes/wvgeometrydoublyperiodic/maskforconjugatefouriercoefficients.html) a mask indicate the components that are redundant conjugates
  + [`maskForNyquistModes`](/classes/wvgeometrydoublyperiodic/maskfornyquistmodes.html) returns a mask with locations of modes that are not fully resolved
+ Utility function
  + [`degreesOfFreedomForComplexMatrix`](/classes/wvgeometrydoublyperiodic/degreesoffreedomforcomplexmatrix.html) a matrix with the number of degrees-of-freedom at each entry
  + [`degreesOfFreedomForRealMatrix`](/classes/wvgeometrydoublyperiodic/degreesoffreedomforrealmatrix.html) a matrix with the number of degrees-of-freedom at each entry
  + [`indicesOfFourierConjugates`](/classes/wvgeometrydoublyperiodic/indicesoffourierconjugates.html) a matrix of linear indices of the conjugate
  + [`isHermitian`](/classes/wvgeometrydoublyperiodic/ishermitian.html) Check if the matrix is Hermitian. Report errors.
  + [`setConjugateToUnity`](/classes/wvgeometrydoublyperiodic/setconjugatetounity.html) set the conjugate of the wavenumber (iK,iL) to 1


---